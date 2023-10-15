#include <Eigen/Eigen>
#include <ImfRgbaFile.h>
#include <iostream>
#include <tbb/tbb.h>
#include <vector>
#include <fstream>
#include "Scene.h"
#include <random>
#include <ImfRgba.h>
#include <ImfArray.h>
#include <cmath>
#include <vector>
#include <shared_mutex>

#define SPHERICAL_ENVMAP
//#include <openexr.h>

//using namespace std;
using namespace Eigen;
#define WIDTH 1366
#define HEIGHT 1024
#define PI 3.1415926f

int gCubeLightWidth = 0;
int gCubeLightHeight = 0;
using namespace std;

template<typename T>
int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

vector<Vector3f> gImage(WIDTH* HEIGHT, Vector3f(0.0f, 0.0f, 0.0f));
shared_ptr< Imf::Array2D<Imf::Rgba>> lightImage;
void RenderGPU();

Scene gScene;


Vector3f gPos(0.f,0.f,0.f);
//Vector3f gLookAt(0.f,0.f,0.f);
Matrix3f gView;

struct BSSRDFTable
{
	const int nRhoSamples, nRadiusSamples;
	unique_ptr<float[]> rhoSamples, radiusSamples;
	unique_ptr<float[]> profile;
	unique_ptr<float[]> rhoEff;
	unique_ptr<float[]> profileCDF;
	BSSRDFTable(int nRhoSamples, int nRadiusSamples)
		: nRhoSamples(nRhoSamples),
		nRadiusSamples(nRadiusSamples),
		rhoSamples(new float[nRhoSamples]),
		radiusSamples(new float[nRadiusSamples]),
		profile(new float[nRadiusSamples * nRhoSamples]),
		rhoEff(new float[nRhoSamples]),
		profileCDF(new float[nRadiusSamples * nRhoSamples]) {}
	inline float EvalProfile(int rhoIndex, int radiusIndex) const {
		return profile[rhoIndex * nRadiusSamples + radiusIndex];
	}
};

unique_ptr<BSSRDFTable> gTable;

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

template<typename Predicate>
int FindInterval(int size, const Predicate& pred)
{
	int first = 0, len = size;
	while (len >0)
	{
		int half = len >> 1, middle = first + half;
		if (pred(middle)) {
			first = middle + 1;
			len -= half + 1;
		}
		else
			len = half;
	}
	return Clamp(first - 1, 0, size - 2);
}
struct Distribution1D {
	Distribution1D(const float* f, int n) : func(f, f + n), cdf(n + 1)
	{
		cdf[0] = 0;
		for (int i = 1; i < n + 1; i++) cdf[i] = cdf[i - 1] + func[i - 1] / n;
		funcInt = cdf[n];
		if (funcInt == 0) {
			for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
		}
		else {
			for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
		}
	}
	int Count() const { return (int)func.size(); }
	float SampleContinuous(float u, float* pdf, int* off = nullptr) const {
		int offset = FindInterval((int)cdf.size(),
			[&](int index) { return cdf[index] <= u; });
		if (off) *off = offset;
		float du = u - cdf[offset];
		if ((cdf[offset + 1] - cdf[offset]) > 0)
		{
			du /= (cdf[offset + 1] - cdf[offset]);
		}
		if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;
		return (offset + du) / Count();
	}
	int SampleDiscrete(float u, float* pdf = nullptr,
		float* uRemapped = nullptr) const {
		int offset = FindInterval((int)cdf.size(),
			[&](int index) { return cdf[index] <= u; });
		if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
		if (uRemapped)
			*uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
	//	if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
		return offset;
	}
	vector<float> func, cdf;
	float funcInt;
};

class Distribution2D {
public:
	// Distribution2D Public Methods
	Distribution2D(const float* data, int nu, int nv);
	Vector2f SampleContinuous(const Vector2f& u, float* pdf) const {
		float pdfs[2];
		int v;
		float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
		float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
		*pdf = pdfs[0] * pdfs[1];
		return Vector2f(d0, d1);
	}
	float Pdf(const Vector2f& p) const {
		int iu = Clamp(int(p[0] * pConditionalV[0]->Count()), 0,
			pConditionalV[0]->Count() - 1);
		int iv =
			Clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
		return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
	}

private:
	// Distribution2D Private Data
	std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
	std::unique_ptr<Distribution1D> pMarginal;
};

Distribution2D::Distribution2D(const float* func, int nu, int nv) {
	pConditionalV.reserve(nv);
	for (int v = 0; v < nv; ++v) {
		// Compute conditional sampling distribution for $\tilde{v}$
		pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
	}
	// Compute marginal sampling distribution $p[\tilde{v}]$
	std::vector<float> marginalFunc;
	marginalFunc.reserve(nv);
	for (int v = 0; v < nv; ++v)
		marginalFunc.push_back(pConditionalV[v]->funcInt);
	pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
}


unique_ptr<Distribution2D> gLightDistribution;

float FresnelMoment1(float eta)
{
	float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta, eta5 = eta4 * eta;
	if (eta < 1)
		return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 + 2.49277f * eta4 - 0.68441f * eta5;
	else
		return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 - 1.27198f * eta4 + 0.12746f * eta5;
}
float FresnelMoment2(float eta)
{
	float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta, eta5 = eta4 * eta;
	if (eta < 1)
		return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 + 0.07883f * eta4 + 0.4860f * eta5;
	else
	{
		float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
		return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 + 458.843f * r_eta + 404.557f * eta - 189.519f * eta2 + 54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
	}
}
/// <summary>
/// 颜色到光学参数（sigma_s, sigma_a)
/// </summary>
/// <param name="fragCoord"></param>
/// <returns></returns>
float BeamDiffusionMS(float sigma_s, float sigma_a, float g, float eta, float r)
{
	float sigmap_s = sigma_s * (1 - g);
	float sigmap_t = sigma_a + sigmap_s;
	float rhop = sigmap_s / sigmap_t;
	float D_g = (2 * sigma_a + sigmap_s) / (3 * sigmap_t * sigmap_t);
	float sigma_tr = sqrt(sigma_a / D_g);
	float fm1 = FresnelMoment1(eta), fm2 = FresnelMoment2(eta);
	float ze = -2 * D_g * (1 + 3 * fm2) / (1 - 2 * fm1);
	float cPhi = .25f * (1 - 2 * fm1), cE = .5f * (1 - 3 * fm2);
	const int nSamples = 100;
	float Ed = 0;
	for (int i=0;i<nSamples;i++)
	{
		float zr = -log(1 - (i + .5f) / nSamples) / sigmap_t;
		float zv = -zr + 2 * ze;
		float dr = sqrt(r * r + zr * zr), dv = sqrt(r * r + zv * zv);
		float phiD = 1.0f / (4.f * PI) / D_g * (exp(-sigma_tr * dr) / dr - exp(-sigma_tr * dv) / dv);
		float EDn = 1.0f / (4.f * PI) * (zr * (1 + sigma_tr * dr) * exp(-sigma_tr * dr) / (dr * dr * dr) - zv * (1 + sigma_tr * dv) * exp(-sigma_tr * dv) / (dv * dv * dv));
		float E = phiD * cPhi + EDn * cE;
		float kappa = 1 - exp(-2 * sigmap_t * (dr + zr));
		Ed += kappa * rhop * rhop * E;
	}
	return Ed / nSamples;
}
float FrDielectric(float cosThetaI, float etaI, float etaT) {
	cosThetaI = min(max(cosThetaI, -1), 1);
	// Potentially swap indices of refraction
	bool entering = cosThetaI > 0.f;
	if (!entering) {
		std::swap(etaI, etaT);
		cosThetaI = std::abs(cosThetaI);
	}

	// Compute _cosThetaT_ using Snell's law
	float sinThetaI = std::sqrt(max((float)0, 1 - cosThetaI * cosThetaI));
	float sinThetaT = etaI / etaT * sinThetaI;

	// Handle total internal reflection
	if (sinThetaT >= 1) return 1;
	float cosThetaT = std::sqrt(max((float)0, 1 - sinThetaT * sinThetaT));
	float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
		((etaT * cosThetaI) + (etaI * cosThetaT));
	float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
		((etaI * cosThetaI) + (etaT * cosThetaT));
	return (Rparl * Rparl + Rperp * Rperp) / 2;
}
inline float PhaseHG(float cosTheta, float g) {
	float denom = 1 + g * g + 2 * g * cosTheta;
	return PI/4.0f * (1 - g * g) / (denom * std::sqrt(denom));
}
float BeamDiffusionSS(float sigma_s, float sigma_a, float g, float eta,
	float r) {
	// Compute material parameters and minimum $t$ below the critical angle
	float sigma_t = sigma_a + sigma_s, rho = sigma_s / sigma_t;
	float tCrit = r * std::sqrt(eta * eta - 1);
	float Ess = 0;
	const int nSamples = 100;
	for (int i = 0; i < nSamples; ++i) {
		// Evaluate single scattering integrand and add to _Ess_
		float ti = tCrit - std::log(1 - (i + .5f) / nSamples) / sigma_t;

		// Determine length $d$ of connecting segment and $\cos\theta_\roman{o}$
		float d = std::sqrt(r * r + ti * ti);
		float cosThetaO = ti / d;

		// Add contribution of single scattering at depth $t$
		Ess += rho * std::exp(-sigma_t * (d + tCrit)) / (d * d) *
			PhaseHG(cosThetaO, g) * (1 - FrDielectric(-cosThetaO, 1, eta)) *
			std::abs(cosThetaO);
	}
	return Ess / nSamples;
}

float gEta = 1.5f;

Vector3f gSigma_t(1000,1000,1000.f);


// float Integrate(float rho)
// {
// 	float r[64];
// 	r[0] = 0;
// 	r[1] = 2.5e-3f;
// 	for (int i=2;i<64;i++)
// 	{
// 		r[i] = r[i - 1] * 1.2f;
// 	}
// 	float values[64];
// 	for (int i=0;i<64;i++)
// 	{
// 		//values[i] = (1.f - exp(-8.0f * i / 63.0f)) / (1.f - exp(-8.f));
// 		values[i]=2.0f*PI*r[i] * BeamDiffusionMS(rho, 1 - rho, 0, gEta, r[i]);
// 	}
// 	float sum = 0.f;
// 	for (int i=0;i<63;i++)
// 	{
// 		sum += (values[i] + values[i + 1]) / 2.0f * (r[i + 1] - r[i]);
// 	}
// 	return sum;
// }

bool CatmullRomWeights(int size, const float* nodes, float x, int* offset,
	float* weights) {
	if (!(x >= nodes[0] && x <= nodes[size - 1])) return false;
	int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
	*offset = idx - 1;
	float x0 = nodes[idx], x1 = nodes[idx + 1];
	float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;
	weights[1] = 2 * t3 - 3 * t2 + 1;
	weights[2] = -2 * t3 + 3 * t2;
	if (idx > 0) {
		float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
		weights[0] = -w0;
		weights[2] += w0;
	}
	else {
		float w0 = t3 - 2 * t2 + t;
		weights[0] = 0;
		weights[1] -= w0;
		weights[2] += w0;
	}
	if (idx + 2 < size) {
		float w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
		weights[1] -= w3;
		weights[3] = w3;
	}
	else {
		float w3 = t3 - t2;
		weights[1] -= w3;
		weights[2] += w3;
		weights[3] = 0;
	}
	return true;
}
Vector3f Sp(float r, const Vector3f& rho)
{
// 	Vector3f Sr;
// 	for (int ch=0;ch<3;ch++)
// 	{
// 		float rOptical = r * gSigma_t[ch];
// 		Sr[ch] = BeamDiffusionMS(rho[ch], 1 - rho[ch], 0, gEta, rOptical);
// 	}
// 	Sr.array() *= gSigma_t.array() * gSigma_t.array();
// 
	Vector3f Sr=Vector3f::Constant(0.f);
	for (int ch = 0; ch < 3; ++ch) {
		// Convert $r$ into unitless optical radius $r_{\roman{optical}}$
		float rOptical = r * gSigma_t[ch];

		// Compute spline weights to interpolate BSSRDF on channel _ch_
		int rhoOffset, radiusOffset;
		float rhoWeights[4], radiusWeights[4];
		if (!CatmullRomWeights(gTable->nRhoSamples, gTable->rhoSamples.get(),
			rho[ch], &rhoOffset, rhoWeights) ||
			!CatmullRomWeights(gTable->nRadiusSamples, gTable->radiusSamples.get(),
				rOptical, &radiusOffset, radiusWeights))
			continue;

		// Set BSSRDF value _Sr[ch]_ using tensor spline interpolation
		float sr = 0;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				float weight = rhoWeights[i] * radiusWeights[j];
				if (weight != 0)
					sr += weight * gTable->EvalProfile(rhoOffset + i, radiusOffset + j);
			}
		}

		// Cancel marginal PDF factor from tabulated BSSRDF profile
		if (rOptical != 0) sr /= 2 * PI * rOptical;
		Sr[ch] = sr;
	}
	// Transform BSSRDF value into world space units
	Sr.array() *= gSigma_t.array() * gSigma_t.array();
	return Sr;
}



float Pdf_Sr(int ch, float r, float rho)
{
	float rOptical = r * gSigma_t[ch];

	int rhoOffset, radiusOffset;
	float rhoWeights[4], radiusWeights[4];
	if (!CatmullRomWeights(gTable->nRhoSamples, gTable->rhoSamples.get(), rho,
		&rhoOffset, rhoWeights) ||
		!CatmullRomWeights(gTable->nRadiusSamples, gTable->radiusSamples.get(),
			rOptical, &radiusOffset, radiusWeights))
		return 0.f;
	float sr = 0, rhoEff = 0;
	for (int i = 0; i < 4; ++i) {
		if (rhoWeights[i] == 0) continue;
		rhoEff += gTable->rhoEff[rhoOffset + i] * rhoWeights[i];
		for (int j = 0; j < 4; ++j) {
			if (radiusWeights[j] == 0) continue;
			sr += gTable->EvalProfile(rhoOffset + i, radiusOffset + j) *
				rhoWeights[i] * radiusWeights[j];
		}
	}
	if (rOptical != 0) sr /= 2 * PI * rOptical;
	return std::max<float>((float)0, sr * gSigma_t[ch] * gSigma_t[ch] / rhoEff);
}

float pdf_Sp(const Vector3f& ss, const Vector3f& ts, const Vector3f& ns, const Intersection& po, const Intersection& pi, const Vector3f& rho)
{
	Vector3f d = po.coords-pi.coords;
	Vector3f dLocal(ss.dot(d), ts.dot(d), ns.dot(d));
	Vector3f nLocal(ss.dot(pi.normal), ts.dot(pi.normal), ns.dot(pi.normal));
	float rProj[3] = { std::sqrt(dLocal.y() * dLocal.y() + dLocal.z() * dLocal.z()),
					  std::sqrt(dLocal.z() * dLocal.z() + dLocal.x() * dLocal.x()),
					  std::sqrt(dLocal.x() * dLocal.x() + dLocal.y() * dLocal.y()) };
	float pdf = 0, axisProb[3] = { .25f, .25f, .5f };
	float chProb = 1.0f / 3;
	for (int axis = 0; axis < 3; ++axis)
		for (int ch = 0; ch < 3; ++ch)
			pdf += Pdf_Sr(ch, rProj[axis], rho[ch]) * std::abs(nLocal[axis]) * chProb * axisProb[axis];
	return pdf;
}






template<typename T>
T lerp(T a, T b, float c)
{
	return a * (1 - c) + b * c;
}

Imf::Rgba operator*(Imf::Rgba a, float b)
{
	Imf::Rgba out;
	out.r = a.r * b;
	out.g = a.g * b;
	out.b = a.b * b;
	out.a = a.a * b;
	return out;
}
Imf::Rgba operator+(Imf::Rgba a, Imf::Rgba b)
{
	Imf::Rgba out;
	out.r = a.r + b.r;
	out.g = a.g + b.g;
	out.b = a.b + b.b;
	out.a = a.a + b.a;
	return out;
}
Vector3f SampleCubeLight(Vector2f uv)
{
	Vector2f coord = uv.array() * Vector2f(gCubeLightWidth - 1, gCubeLightHeight - 1).array();
	Vector2i coordInt = coord.cast<int>();
	Vector2f fract = coord - coordInt.cast<float>();
	Imf::Rgba sample1 = (*lightImage)[coordInt.y()][coordInt.x()];
	Imf::Rgba sample2 = (*lightImage)[coordInt.y()][min(coordInt.x() + 1, gCubeLightWidth - 1)];
	Imf::Rgba sample3 = (*lightImage)[min(coordInt.y() + 1, gCubeLightHeight - 1)][coordInt.x()];
	Imf::Rgba sample4 = (*lightImage)[min(coordInt.y() + 1, gCubeLightHeight - 1)][min(coordInt.x() + 1, gCubeLightWidth - 1)];
	Imf::Rgba color = lerp(lerp(sample1, sample2, fract.x()), lerp(sample3, sample4, fract.x()), fract.y());
	//Imf::Rgba color = (*lightImage)[(int)(uv.y() * (gCubeLightHeight-1))][(int)(uv.x() * (gCubeLightWidth-1))];
	return Vector3f(color.r, color.g, color.b);
}

Vector2f SphericalLightToUV(Vector3f p, float *pOutSinTheta=nullptr)
{
	Vector2f uv;

	Matrix4f mat1, mat2;
	mat1 << 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1;
	mat2 << -0.224951, 0.0, -0.97437, 0.0, -0.97437, 0.0, 0.224951, 0.0, 0.0, 1.0, 0.0, 8.87, 0.0, 0.0, 0.0, 1.0;
	Vector4f pp = mat1 * mat2 * Vector4f(p.x(), p.y(), p.z(), 0.0);
	p.x() = pp.x();
	p.y() = pp.y();
	p.z() = pp.z();
	float phi = fmod((atan2(p.y(), p.x()) + 2 * PI), 2.0f * PI);
	float theta = acos(p.z());
	uv.array() = Vector2f(phi, theta).array() / Vector2f(2. * PI, PI).array();
	if (pOutSinTheta)
		*pOutSinTheta = sin(theta);
	return uv;
}
Vector3f SampleCubeLight(Vector3f p)
{
	// transform p 

	Vector2f uv;
#ifdef SPHERICAL_ENVMAP
	uv = SphericalLightToUV(p);
#else
	float total = abs(p.x()) + abs(p.y()) + abs(p.z());
	p /= total;

	if (p.y() > 0.)
		uv = Vector2f(p.x(), p.z()) * 0.5 + Vector2f(0.5, 0.5);
	else
		//uv = Vector2f(p.x(), p.z()) * 0.5 + Vector2f(0.5, 0.5);
		uv.array() = (Vector2f(1.,1.) - Vector2f(abs(p.z()), abs(p.x()))).array() * Vector2f(sign(p.x()), sign(p.z())).array() * 0.5 + Vector2f(0.5, 0.5).array();

#endif
	//uv.array() = Vector2f(p.x(), p.y()).array() / Vector2f(WIDTH, HEIGHT).array();
	return SampleCubeLight(uv);
	//fragColor = vec4(texture(iChannel0, texc));
}
Vector3f CosineSampleHemisphere(float& pdf)
{
	Vector2f _rand(get_random_float(), get_random_float());
	//Vector2f _rand(0.0, 1.0);
	float Phi = 2 * PI * _rand.x();
	float CosTheta = sqrt(_rand.y());
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	Vector3f H;
	H.x() = SinTheta * cos(Phi);
	H.y() = SinTheta * sin(Phi);
	H.z() = CosTheta;

	pdf = CosTheta * (1.0 / PI);

	return H;
}
#define DegreeToRadian(d)  (d/180.f*PI)

Vector3f gSigmaS(0.74, 0.88, 1.01);
Vector3f gSigmaA(0.032, 0.17, 0.48);


// float get_random_float1()
// {
// 	static float s_inc = 0.0f;
// 	s_inc += 0.001f;
// 	if (s_inc > 1.0f)
// 		s_inc = 0.0f;
// 	return s_inc;
// 
// }

float lds(uint32_t i)
{
	const uint32_t base = 2;
	const float invbase = 1.0f / (float)base;
	uint32_t reversedDigits = 0;
	float invBaseN = 1;
	while (i)
	{
		uint32_t next = i / base;
		uint32_t digit = i - next * base;
		reversedDigits = reversedDigits * base + digit;
		invBaseN *= invbase;
		i = next;
	}
	return reversedDigits * invBaseN;
}

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

#define SPP 4

float Sw(float cosTheta)
{
	float c = 1 - 2 * FresnelMoment1(1 / gEta);
	return (1 - FrDielectric(cosTheta, 1, gEta)) / (c * PI);
}

std::atomic_int gSwCount = 0;

Vector3f SampleOneLight(const Intersection& inter)
{
	Vector3f Color = Vector3f::Constant(0);
	{
		Vector2f uv(get_random_float(), get_random_float());
		float pdf;
		Vector2f sample = gLightDistribution->SampleContinuous(uv, &pdf);
		float theta = sample[1] * PI, phi = sample[0] * 2 * PI;
		float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
		float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
		Matrix4f mat1, mat2;
		mat1 << 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1;
		mat2 << -0.224951, 0.0, -0.97437, 0.0, -0.97437, 0.0, 0.224951, 0.0, 0.0, 1.0, 0.0, 8.87, 0.0, 0.0, 0.0, 1.0;
		Matrix4f matInverse = (mat1 * mat2).inverse();
		Vector4f _dir = matInverse * Vector4f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta, 0.f);
		Vector3f dir(_dir.x(), _dir.y(), _dir.z());
		pdf /= (2 * PI * PI * sinTheta);
		if (sinTheta == 0) pdf = 0;
		if (pdf > 0)
		{
			Ray _ray(inter.coords + inter.normal * 0.00001f, dir);
			Intersection inter3 = gScene.intersect(_ray);
			if (!inter3.happened)
			{
				Vector3f li = SampleCubeLight(sample);
				float ndotl = max(0.f, dir.dot(inter.normal));
				// include sw for bssrdf
				//Vector3f bsdf = inter.color.array() / PI;
				Vector3f bsdf = Vector3f::Constant(Sw(ndotl)) * gEta * gEta;
				
				float brdfpdf = ndotl / PI;
				float weight = PowerHeuristic(1, pdf, 1, brdfpdf);
				
				Color.array() += bsdf.array() * li.array() * ndotl / pdf * weight;
	
			}
		}
	}



	{
		float pdf;
		Vector3f dir = CosineSampleHemisphere(pdf);
		Vector3f ldir = dir;
		Vector3f bsdf = Vector3f::Constant(Sw(dir.z())) * gEta * gEta;
		//Vector3f bsdf = inter.color.array() / PI;

		//transform to world
		Vector3f ss(1,0,0);
		Vector3f ts = inter.normal.cross(ss);
		if (ts.norm() == 0)
		{
			ss = Vector3f(0, 1, 0);
			ts = inter.normal.cross(ss);
		}
		ts.normalize();
		ss = ts.cross(inter.normal).normalized();

		dir = ss * dir.x() + ts * dir.y() + inter.normal * dir.z();

		Ray _ray(inter.coords + inter.normal * 0.00001f, dir);
		Intersection inter3 = gScene.intersect(_ray);
		if (!inter3.happened)
		{
			Vector3f li = SampleCubeLight(dir);
			float SinTheta = 0.f;
			Vector2f uv = SphericalLightToUV(dir, &SinTheta);

			float lightpdf = gLightDistribution->Pdf(uv) / (2 * PI * PI * SinTheta);
			float weight = PowerHeuristic(1, pdf, 1, lightpdf);
			float ndotl = max(0.f, dir.dot(inter.normal));
			Color.array() += /*inter.color.array() **/ li.array() * bsdf.array() * ndotl /pdf * weight;
			gSwCount++;
		}
	}
	return Color;
}

inline Vector3f Reflect(const Vector3f& wi, const Vector3f& n)
{
	return wi - wi.dot(n) * n * 2;
}
inline Vector3f Faceforward(const Vector3f& n, const Vector3f& v) {
	return (n.dot(v) < 0.f) ? -n : n;
}
bool Refract(const Vector3f& wi, const Vector3f& n, float eta, Vector3f& wt)
{
	float cosi = wi.dot(n);
	float sini2 = 1 - cosi * cosi;
	float sino2 = sini2 * eta*eta;
	if (sino2 >= 1) return false;
	Vector3f wip = cosi * n;
	Vector3f wih = wi - wip;
	Vector3f woh = -wih * eta;
	float coso = sqrt(1 - sino2);
	Vector3f wop = -n * coso;
	wt = (wop + woh).normalized();
	return true;
}

Vector3f SampleSpecular(const Vector3f &wo, const Vector3f &normal, Vector3f& dir, float& pdf, bool &bTransmission)
{
	float u = get_random_float();
	float F = FrDielectric(wo.dot(normal), 1, gEta);
	if (u < F) {
		// Compute specular reflection for _FresnelSpecular_

		// Compute perfect specular reflection direction
		dir = Reflect(-wo,normal);
		bTransmission = false;
		pdf = F;
		return Vector3f::Constant(F) / max(0, dir.dot(normal));
	}
	else {
		// Compute specular transmission for _FresnelSpecular_

		// Figure out which $\eta$ is incident and which is transmitted
		bool entering = wo.dot(normal) > 0;
		float etaI = entering ? 1 : gEta;
		float etaT = entering ? gEta : 1;

		// Compute ray direction for specular transmission
		if (!Refract(wo, Faceforward(normal, wo), etaI / etaT, dir))
			return Vector3f::Constant(0);
		float ft = (1 - F);

		// Account for non-symmetry with transmission to different medium
//		if (mode == TransportMode::Radiance)
			ft *= (etaI * etaI) / (etaT * etaT);
		bTransmission = true;
		pdf = 1 - F;
		//ft *= 0.5f;
		return Vector3f::Constant(ft) / abs(dir.dot(normal));
	}
}



float SampleCatmullRom2D(int size1, int size2, const float* nodes1,
	const float* nodes2, const float* values,
	const float* cdf, float alpha, float u, float* fval = nullptr,
	float* pdf = nullptr) {
	int offset;
	float weights[4];
	if (!CatmullRomWeights(size1, nodes1, alpha, &offset, weights)) return 0;
	auto interpolate = [&](const float* array, int idx) {
		float value = 0;
		for (int i = 0; i < 4; ++i)
			if (weights[i] != 0)
				value += array[(offset + i) * size2 + idx] * weights[i];
		return value;
		};
	float maximum = interpolate(cdf, size2 - 1);
	u *= maximum;
	int idx =
		FindInterval(size2, [&](int i) { return interpolate(cdf, i) <= u; });
	
	float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
	float x0 = nodes2[idx], x1 = nodes2[idx + 1];
	float width = x1 - x0;
	float d0, d1;
	u = (u - interpolate(cdf, idx)) / width;
	if (idx > 0)
		d0 = width * (f1 - interpolate(values, idx - 1)) / (x1 - nodes2[idx - 1]);
	else
		d0 = f1 - f0;
	if (idx + 2 < size2)
		d1 = width * (interpolate(values, idx + 2) - f0) /(nodes2[idx + 2] - x0);
	else
		d1 = f1 - f0;
	float t;
	if (f0 != f1)
		t = (f0 - std::sqrt(std::max<float>(0, f0 * f0 + 2 * u * (f1 - f0)))) / (f0 - f1);
	else
		t = u / f0;

	float a = 0, b = 1, Fhat, fhat;
	while (true)
	{
		if (!(t >= a && t <= b)) t = 0.5f * (a + b);
		Fhat = t * (f0 +
			t * (.5f * d0 +
				t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
					t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
		fhat = f0 +
			t * (d0 +
				t * (-2 * d0 - d1 + 3 * (f1 - f0) +
					t * (d0 + d1 + 2 * (f0 - f1))));

		if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;
		if (Fhat - u < 0)
			a = t;
		else
			b = t;
		t -= (Fhat - u) / fhat;
	}
	if (fval) *fval = fhat;
	if (pdf) *pdf = fhat / maximum;
	return x0 + t * width;
}

float Sample_Sr(int ch, float u, float rho) {
	if (gSigma_t[ch] == 0) return -1;
	return SampleCatmullRom2D(gTable->nRhoSamples, gTable->nRadiusSamples,
		gTable->rhoSamples.get(), gTable->radiusSamples.get(),
		gTable->profile.get(), gTable->profileCDF.get(),
		rho, u) / gSigma_t[ch];
}

Vector3f Scatter(const Intersection& inter, Intersection& inter1, float& pdf)
{


	Vector3f vx, vy, vz;
	float u = get_random_float();
	//float u = 0.5;
	Vector3f ss = inter.normal + Vector3f(0,0,0.1f);
	ss.normalize();
	Vector3f ts = inter.normal.cross(ss).normalized();
	ss = ts.cross(inter.normal).normalized();
	if (u<.5f)
	{
		vx = ss;
		vy = ts;
		vz = inter.normal;
		u *= 2;

	}
	else if (u<0.75f)
	{
		vx = ts;
		vy = inter.normal;
		vz = ss;
		u = (u - .5f) * 4;

	}
	else
	{
		vx = inter.normal;
		vy = ss;
		vz = ts;
		u = (u - .75f) * 4;

	}

	int ch = min(2, (int)(u * 3));
	u = u * 3 - ch;
	gSigma_t = gSigmaS + gSigmaA;
	float sigmat = gSigma_t[ch];
	Vector3f rhos = gSigmaS.array() / gSigma_t.array();
	float rho = gSigmaS[ch] / sigmat;
	float r = Sample_Sr(ch, get_random_float(), rho);
	if (r < 0) return Vector3f::Constant(0.f);
	float phi = 2 * PI * get_random_float();


	float r_max = Sample_Sr(ch, 0.999f, rho);
	if (r >= r_max) return Vector3f::Constant(0.f);

	//sample
	float l = 2 * sqrt(r_max * r_max - r * r);


	Ray _ray(inter.coords + r * (vx * cos(phi) + vy * sin(phi)) - l * vz * 0.5f, vz);
	vector<Intersection> intersects;
	Intersection _inter;
	float _l = 0;
	Vector3f origin = _ray.origin;
	while ((_inter = gScene.intersect(_ray)).happened && _l < l)
	{
		intersects.push_back(_inter);
		_l = (_inter.coords - origin).norm();
		_ray.origin = _inter.coords + _ray.direction * 0.00001;
	}
	int size = intersects.size();
	if (size != 0)
	{
		int selected = min(u * size, size - 1);
		inter1 = intersects[selected];
		pdf = pdf_Sp(ss,ts,inter.normal, inter, inter1, rhos) / size;
		return Sp((inter.coords - inter1.coords).norm(), rhos);

		//inter2.coords
		//float pdf;
		//Vector3f dir = CosineSampleHemisphere(pdf);
		////transform to world
		//Vector3f ss = inter2.normal + Vector3f(0, 0, 0.1f);
		//ss.normalize();
		//Vector3f ts = inter2.normal.cross(ss).normalized();
		//ss = ts.cross(inter2.normal).normalized();
		//dir = ts * dir.x() + inter2.normal * dir.y() + ss * dir.z();


		//float invpdf = gRadicalInt * size / sigmat/sigmat;
		//Vector3f ptoi = inter1.coords - inter2.coords;

		//float scale = Fro * invpdf * (1 - FrDielectric(inter2.normal.dot(dir), 1.0, gEta))/(gC*PI);

		//float dist = ptoi.norm();
		//Ray _ray(inter2.coords + inter2.normal * 0.00001f, dir);
		//Intersection inter3 = gScene.intersect(_ray);
		//if (!inter3.happened)
		//{
		//	//Vector3f brdf = inter2.color.array() / 3.1415926f;
		//	Vector3f le = SampleCubeLight(dir);
		//	Color.array() += 100.0f * scale * inter2.color.array() * le.array();
		//	//Color.array() += 5000.0f * scale * brdf.array() * max(0.f, -li.dot(inter1.normal)) * max(0.f, li.dot(inter2.normal)) / (dist * dist) * inter3.emit.array() / pdf;
		//}
	}
	return Vector3f::Constant(0);
}
std::atomic_int gZeroCount = 0;
std::atomic_int gPassedCount = 0;

Vector3f Render(Vector2f fragCoord)
{
	Vector2f uv = fragCoord.array() / Vector2f(WIDTH, HEIGHT).array() * 2.0f - 1.0f;
	uv.y() = -uv.y();
	
	uv *= tan(DegreeToRadian(28.8415038750464 / 2.f));
	//uv *= tan(PI / 3.f);

	uv.x() *= (float)WIDTH / HEIGHT;
	std::mutex color_visit;
	Vector3f dir = gView * Vector3f(uv.x(), uv.y(), -1.0f).normalized();
	//Vector3f dir = Vector3f(uv.x(), uv.y(), 1.0f).normalized();
	Ray r(gPos, dir);
	Intersection inter = gScene.intersect(r);
	if (inter.happened)
	//if (false)
	{
		if (inter.bEmit)
		{
			return inter.emit;
		}
		else
		{
			Vector3f Color(0.f, 0.f, 0.f);// = inter.color;
			vector<Light*> lights = gScene.GetAllLights();
			
		
			//VLOG(2) << "Updated beta = " << beta;

			//float Fro = 1 - FrDielectric(-r.direction.dot(inter.normal),1.0f,gEta);

			const int spp = 1; // SPP* SPP;
			tbb::parallel_for(tbb::blocked_range<int>(0, spp), [&](tbb::blocked_range<int> i) {
				for (int j = i.begin(); j != i.end(); j++)
					//for(int j=0;j<spp;++j)
				{
					float pdf;
					//Vector3f dir = CosineSampleHemisphere(pdf);

					Vector3f wi;
					bool bTransmission = false;
					Vector3f f = SampleSpecular(-dir, inter.normal, wi, pdf, bTransmission);

// 					int i = 0;
// 					if (fragCoord.x() == 600 && fragCoord.y() == 400)
// 					{
// 						i++;
// 					}

					//if (f.IsBlack() || pdf == 0.f) break;
					Vector3f beta = f * abs(wi.dot(inter.normal)) / pdf;
					//beta = 1.0f;
					// beta FrDeletric (1-Fr(theta)  for bssrdf
					//beta = Vector3f::Constant(1.0f);
					if (bTransmission)
					{
						//transmission and bssrdf
						Intersection inter1;
						pdf = 0;
						Vector3f S = Scatter(inter, inter1, pdf);
						if (S.norm() == 0.f || pdf == 0)
						{
							//beta.array() *= Vector3f::Constant(0.0).array();
							gZeroCount++;
						}
						else
						{
								//SampleOneLight(inter, color_visit, Color);
							beta.array() *= S.array() / pdf;
							gPassedCount++;
							std::lock_guard<std::mutex> lg(color_visit);
							Color.array() += beta.array() * SampleOneLight(inter1).array();
						}
						
					}
					else
					{
						// reflection
						std::lock_guard<std::mutex> lg(color_visit);
						Color.array() += beta.array() * SampleCubeLight(wi).array();
					}

					//for (Light* light : lights)
					//{
					//	Intersection inter1;
						//float pdf;
						//light->Sample(inter1, pdf);

			
						// mis 
						
					
						//spherical sample
						//Vector2f uv(float(j% SPP)/SPP,float(j/SPP)/SPP);
						
						//Vector2f uv(lds(nCount),(float)nCount/(WIDTH*HEIGHT));







					/*	Vector3f ptoi = inter1.coords - inter.coords;
						float dist = ptoi.norm();
						Vector3f li = ptoi.normalized();
						Ray _ray(inter.coords + inter.normal * 0.00001f, li);
						Intersection inter3 = gScene.intersect(_ray);
						float len = (inter1.coords - inter3.coords).norm();
						if (inter3.happened && len < 1e-3 )
						{
							Vector3f brdf = inter.color.array() / 3.1415926f;

							Color.array() += 5.0f * brdf.array() * max(0.f, -li.dot(inter1.normal)) * max(0.f, li.dot(inter.normal)) / (dist * dist) * inter3.emit.array() / pdf;
						}*/
					//}
				}
			});
			Color /= spp;
// 			for (int i=0;i<3;i++)
// 			{
// 				Color(i) = min(Color(i), 1.f);
// 			}

			//float NdotL = max(0.f,inter.normal.dot(Vector3f(1.f, 1.f, -1.f).normalized()));
			return Color;
		}
	}
	else
		return SampleCubeLight(dir);
}

float IntegrateCatmullRom(int n, const float* x, const float* values,
	float* cdf) {
	float sum = 0;
	cdf[0] = 0;
	for (int i = 0; i < n - 1; ++i) {
		// Look up $x_i$ and function values of spline segment _i_
		float x0 = x[i], x1 = x[i + 1];
		float f0 = values[i], f1 = values[i + 1];
		float width = x1 - x0;

		// Approximate derivatives using finite differences
		float d0, d1;
		if (i > 0)
			d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
		else
			d0 = f1 - f0;
		if (i + 2 < n)
			d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
		else
			d1 = f1 - f0;

		// Keep a running sum and build a cumulative distribution function
		sum += ((d0 - d1) * (1.f / 12.f) + (f0 + f1) * .5f) * width;
		cdf[i + 1] = sum;
	}
	return sum;
}

void InitBSSRDFTable()
{
	gTable = make_unique<BSSRDFTable>(100, 64);
	BSSRDFTable* t = gTable.get();
	t->radiusSamples[0] = 0;
	t->radiusSamples[1] = 2.5e-3f;
	for (int i = 2; i < t->nRadiusSamples; ++i)
		t->radiusSamples[i] = t->radiusSamples[i - 1] * 1.2f;

	// Choose albedo values of the diffusion profile discretization
	for (int i = 0; i < t->nRhoSamples; ++i)
		t->rhoSamples[i] =
		(1 - std::exp(-8 * i / (float)(t->nRhoSamples - 1))) /
		(1 - std::exp(-8));
	tbb::parallel_for(tbb::blocked_range<int>(0,t->nRhoSamples,t->nRhoSamples/16), [&](tbb::blocked_range<int> b) {
		// Compute the diffusion profile for the _i_th albedo sample

		// Compute scattering profile for chosen albedo $\rho$
		for (int i = b.begin(); i<b.end();i++)
		{
			for (int j = 0; j < t->nRadiusSamples; ++j) {
				float rho = t->rhoSamples[i], r = t->radiusSamples[j];
				t->profile[i * t->nRadiusSamples + j] =
					2 * PI * r * (BeamDiffusionSS(rho, 1 - rho, 0.f, gEta, r) 
						/* +BeamDiffusionMS(rho, 1 - rho, 0.f, gEta, r)*/);
			}

			// Compute effective albedo $\rho_{\roman{eff}}$ and CDF for importance
			// sampling
			t->rhoEff[i] =
				IntegrateCatmullRom(t->nRadiusSamples, t->radiusSamples.get(),
					&t->profile[i * t->nRadiusSamples],
					&t->profileCDF[i * t->nRadiusSamples]);
		}
	},tbb::static_partitioner());
}

void main()
{	
	InitBSSRDFTable();
	gSigmaS *= 20;
	gSigmaA *= 20;
 	Imf::RgbaInputFile file("H:\\user\\Desktop\\sssdragon\\textures\\envmap.exr");
 	Imath::Box2i dw = file.dataWindow();
 	gCubeLightWidth = dw.max.x - dw.min.x + 1;
 	gCubeLightHeight = dw.max.y - dw.min.y + 1;
	lightImage = make_shared<Imf::Array2D<Imf::Rgba>>(gCubeLightHeight, gCubeLightWidth);
// 	
// 
 	file.setFrameBuffer(&(*lightImage)[0][0] - dw.min.x - dw.min.y * gCubeLightWidth, 1, gCubeLightWidth);
 	file.readPixels(dw.min.y, dw.max.y);

	vector<float> img;
	img.resize(gCubeLightWidth * gCubeLightHeight);
	const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
	for (int i=0;i<gCubeLightHeight;++i)
	{
		for (int j=0;j<gCubeLightWidth;++j)
		{
			img[i * gCubeLightWidth + j] = YWeight[0] * (*lightImage)[i][j].r + YWeight[1] * (*lightImage)[i][j].g + YWeight[2] * (*lightImage)[i][j].b;
		}
	}
	gLightDistribution.reset(new Distribution2D(&img[0], gCubeLightWidth, gCubeLightHeight));


 
	gScene.LoadScene("H:\\User\\Desktop\\sssdragon\\geometry\\dragon.obj");
	//gScene.LoadScene("D:\\User\\Desktop\\untitled.obj");
//	gScene.AddCubemapLight("");
//	gScene.LoadScene("G:\\desktop\\Assignment7\\models\\bunny\\bunny.obj");

//	auto lookat = gScene.GetObject(0)->getBounds().GetCenter();
	//float radius = gScene.GetObject(0)->getBounds().GetExtent().norm() / 3.0f;
	Vector3f lookat;

//	gScene.Remove(0);
//	gScene.CreateBVH();
	gPos = Vector3f(3.69558, -3.46243, 3.25463);
	lookat = Vector3f(3.04072, - 2.85176, 2.80939);
	Vector3f Up = Vector3f(-0.317366, 0.312466, 0.895346);


	//gPos = lookat - Vector3f(1.2f, -0.7f, 1).normalized() * radius / sin(DegreeToRadian(28.8415038750464 / 2.f));
	//gPos = lookat - Vector3f(0,0,1).normalized() * radius / sin(PI / 4.f);
 	Vector3f Look = (lookat - gPos).normalized();
 	Vector3f Right = Up.cross(Look).normalized();
 	Look = Up.cross(Right).normalized();
 	gView << Right,Up,Look;
	Matrix3f scaleM;
	scaleM << -1, 0, 0, 0, 1, 0, 0, 0, 1;
	gView = gView * scaleM;

	tbb::parallel_for(tbb::blocked_range2d<int>(0, HEIGHT, 0, WIDTH), [&](tbb::blocked_range2d<int> i) {
		{

			for (auto it = i.rows().begin(); it != i.rows().end(); it++)
			{
				for (auto it2 = i.cols().begin(); it2 != i.cols().end(); it2++)
				{
					gImage[it * WIDTH + it2] = Render(Vector2f(it2, it));
				}
			}
		}
	});

	std::cout << "零点数量：" << gZeroCount << std::endl;
	std::cout << "通过点数量：" << gPassedCount << std::endl;
	std::cout << "零点比例：" << (double)gZeroCount/(gZeroCount+gPassedCount) << std::endl;
	std::cout << "SwCount数量：" << gSwCount << std::endl;
	std::cout << "SwCount比例：" << (double)gSwCount / gPassedCount << std::endl;


// 	ofstream ofs("image.ppm");
// 	
// 	ofs << "P3" << endl;
// 	ofs << WIDTH << " "<< HEIGHT << endl;
// 	ofs << 255 << endl;
// 	for (int j=0;j<HEIGHT;j++)
// 	{
// 		for (int i=0;i<WIDTH;i++)
// 		{
// 			int x = min(1.0f,gImage[j * WIDTH + i].x()) * 255.0f;
// 			int y = min(1.0f,gImage[j * WIDTH + i].y()) * 255.0f;
// 			int z = min(1.0f,gImage[j * WIDTH + i].z()) * 255.0f;
// 			ofs <<  x << '\t' <<  y << '\t' << z  << endl;
// 		}
// 	}
// 	ofs.close();

	Imf::Array2D<Imf::Rgba> pixels(HEIGHT,WIDTH);
	for (int j=0;j<HEIGHT;j++)
	{
		for (int i=0;i<WIDTH;i++)
		{
			Imath::half r = gImage[j * WIDTH + i].x();
			Imath::half g = gImage[j * WIDTH + i].y();
			Imath::half b = gImage[j * WIDTH + i].z();
			pixels[j][i] = Imf::Rgba(r, g, b);
		}
	}
	//Imf::Array2D<Imf:Rgba> pixels(HEIGHT, WIDTH);
	Imf::RgbaOutputFile outfile("image.exr", WIDTH, HEIGHT);
	outfile.setFrameBuffer(&pixels[0][0], 1, WIDTH);
	outfile.writePixels(HEIGHT);

// 	Eigen::Vector2d v1;
// 	v1 << 12,23;
// 	Assimp::Importer import;
// 	const char* path = "D:\\User\\Desktop\\Vulkan\\assets\\models\\sponza\\sponza.gltf";
// 	const aiScene * scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);

}