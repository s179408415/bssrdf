#include "mesh.h"


bool Triangle::intersect(const Ray& ray)
{
	return false;
}

bool Triangle::intersect(const Ray&, float&, uint32_t&) const
{
	return false;
}

Intersection Triangle::getIntersection(Ray _ray)
{
 	//Vector3f pvec = _ray.direction.cross(e2);
 	//float det = e1.dot(pvec);
 	//if (fabs(det) < 0.0000001)
 	//	return inter;
 	//float det_inv = 1.f / det;
 	//Vector3f tvec = _ray.origin - v1;
 	//float u = tvec.dot(pvec) * det_inv;
 	//if (u < 0 || u>1)
 	//	return inter;
 	//Vector3f qvec = tvec.cross(e1);
 	//float v = _ray.direction.dot(qvec) * det_inv;
 	//if (v < 0 || u + v>1)
 	//	return inter;
 	//float t_tmp = e2.dot(qvec) * det_inv;
 	//inter.happened = true;
 	//inter.obj = this;
 	//inter.distance = t_tmp;
 
 	//inter.normal = (1 - u - v) * nv1 + u * nv2 + v * nv3;
 	//inter.coords = _ray(t_tmp);
 	//inter.bEmit = hasEmit();
 	//if (hasEmit())
 	//	inter.emit = color;
 	//else
 	//	inter.color = color;


	Intersection inter;
	Vector3f e1 = v2 - v1;
	Vector3f e2 = v3 - v1;
	Vector3f n = e1.cross(e2);
	n.normalize();
// 	float dotval = _ray.direction.dot(n);
// 	if (dotval > 0)
// 		return inter;


	float t_tmp = (v1.dot(n) - _ray.origin.dot(n)) / _ray.direction.dot(n);
	if (t_tmp < 0)
		return inter;
	Vector3f pp = _ray(t_tmp);
	float signArea = (v2 - v1).cross(v3 - v2).dot(n);
	float u = (v3 - v2).cross(pp - v3).dot(n) / signArea;
	float v = (v1 - v3).cross(pp - v1).dot(n) / signArea;
	if (u >= 0 && v >= 0 && u + v <= 1)
	{
		inter.happened = true;
		inter.obj = this;
		inter.distance = t_tmp;

		inter.normal = (u * nv1 + v * nv2 + (1 - u - v) * nv3).normalized();
		inter.coords = pp;
		inter.bEmit = hasEmit();
		if (hasEmit())
			inter.emit = color;
		else
			inter.color = color;

	}

// 	Vector3f _v1 = v1 - _ray.origin;
// 	Vector3f _v2 = v2 - _ray.origin;
// 	Vector3f _v3 = v3 - _ray.origin;
// 
// 	//int maxidx = (int)((float*)std::max_element(&_ray.origin, &_ray.origin + 3)-(float*)(&_ray.origin))/sizeof(float);
// 	vector<float> vals = { _ray.direction.x(), _ray.direction.y(), _ray.direction.z() };
// 	int maxidx = std::max_element(vals.begin(), vals.end())-vals.begin();
// 	int xidx = (maxidx + 1) % 3;
// 	int yidx = (xidx + 1) % 3;
// 	_v1 = Vector3f(_v1[xidx], _v1[yidx], _v1[maxidx]);
// 	_v2 = Vector3f(_v2[xidx], _v2[yidx], _v2[maxidx]);
// 	_v3 = Vector3f(_v3[xidx], _v3[yidx], _v3[maxidx]);
// 	Vector3f rd(_ray.direction[xidx], _ray.direction[yidx], _ray.direction[maxidx]);
// 	float sx = -rd.x() / rd.z();
// 	float sy = -rd.y() / rd.z();
// 	float sz = 1.f / rd.z();
// 	_v1.x() += sx * _v1.z();
// 	_v1.y() += sy * _v1.z();
// 	_v2.x() += sx * _v2.z();
// 	_v2.y() += sy * _v2.z();
// 	_v3.x() += sx * _v3.z();
// 	_v3.y() += sy * _v3.z();
// 
// 	float e0 = _v2.x() * _v3.y() - _v2.y() * _v3.x();
// 	float e1 = _v3.x() * _v1.y() - _v3.y() * _v1.x();
// 	float e2 = _v1.x() * _v2.y() - _v1.y() * _v2.x();
// 	Intersection inter;
// 	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
// 		return inter;
// 	float det = e0 + e1 + e2;
// 	if (det >= 0) return inter;
// 
// 	_v1.z() *= sz;
// 	_v2.z() *= sz;
// 	_v3.z() *= sz;
// 	float tScaled = e0 * _v1.z() + e1 * _v2.z() + e2 * _v3.z();
// 	if (det < 0 && (tScaled >= 0))
// 		return inter;
// 	else if (det > 0 && (tScaled <= 0))
// 		return inter;
// 
// 	float invdet = 1 / det;
// 	float b0 = e0 * invdet;
// 	float b1 = e1 * invdet;
// 	float b2 = e2 * invdet;
// 	float t = tScaled * invdet;
// 
// // 	max(max(_ray.direction.y(), _ray.direction.x()), _ray.direction.z());
// // 	if (u >= 0 && v >= 0 && u + v <= 1)
// // 	{
// 		inter.happened = true;
// 		inter.obj = this;
// 		inter.distance = t;
// 
// 		inter.normal = b0 * nv1 + b1 * nv2 + b2 * nv3;
// 		inter.coords = _ray(t);
// 		inter.bEmit = hasEmit();
// 		if (hasEmit())
// 			inter.emit = color;
// 		else
// 			inter.color = color;
// 
// 	//}
 	return inter;
 }

void Triangle::getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&, Vector2f&) const
{

}

Vector3f Triangle::evalDiffuseColor(const Vector2f&) const
{
	return Vector3f(0.f, 0.f, 0.f);
}

AABB Triangle::getBounds() const
{
	vector<Vector3f> v{ v1, v2, v3 };
	return AABB(v);
}

float Triangle::getArea()
{
	Vector3f e1 = v2 - v1;
	Vector3f e2 = v3 - v1;
	return e1.cross(e2).norm() / 2.0f;
}

void Triangle::Sample(Intersection& pos, float& pdf)
{
	float x = std::sqrt(get_random_float()), y = get_random_float();
	pos.coords = v1 * (1.0f - x) + v2 * (x * (1.0f - y)) + v3 * (x * y);
	Vector3f e1 = v2 - v1;
	Vector3f e2 = v3 - v1;
	pos.normal = e1.cross(e2).normalized();
	pdf = 1.0f / getArea();
}

bool Triangle::hasEmit()
{
	return bEmit;
}

ObjectType Triangle::getType() const
{
	return OT_TRIANGLE;
}

Mesh::Mesh(const vector<Vector3f>& vertices, const vector<Vector3f>& normals, const vector<int>& indices, Vector3f color, bool bEmit)
	//: mAABB(aabb)
{

	int numTris = indices.size() / 3;
	Triangle *triangle = nullptr;
	for (int i=0;i<numTris;i++)
	{
		triangle = new Triangle();
		triangle->v1 = vertices[indices[i * 3]];
		triangle->v2 = vertices[indices[i * 3+1]];
		triangle->v3 = vertices[indices[i * 3+2]];
		triangle->nv1 = normals[indices[i * 3]];
		triangle->nv2 = normals[indices[i * 3+1]];
		triangle->nv3 = normals[indices[i * 3+2]];
		triangle->color = color;
		triangle->bEmit = bEmit;
		mTriangles.push_back(triangle);
	}
	mBVH = BVH::CreateBVH(mTriangles);
}

Mesh::~Mesh()
{
	for (auto &tri : mTriangles)
	{
		delete tri;
	}
	if (mBVH)
		delete mBVH;
}

bool Mesh::intersect(const Ray& ray)
{

	return false;

}

bool Mesh::intersect(const Ray&, float&, uint32_t&) const
{
	return false;
}

Intersection Mesh::getIntersection(Ray _ray)
{
	Intersection intersec;
	if (mBVH)
	{
		intersec = mBVH->GetIntersection(_ray);
	}
	return intersec;
}

void Mesh::getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&, Vector2f&) const
{

}

Vector3f Mesh::evalDiffuseColor(const Vector2f&) const
{
	return Vector3f(0.f, 0.f, 0.f);
}

AABB Mesh::getBounds() const
{
	if (mBVH)
	{
		return mBVH->mAABB;
	}
	return AABB();
}

float Mesh::getArea()
{
	float area = 0.f;
	for (int i =0;i<mTriangles.size();i++)
	{
		area += mTriangles[i]->getArea();
	}
	return area;
}

void Mesh::Sample(Intersection& pos, float& pdf)
{
	float total_area = getArea();
	float acc_area = 0.f;
	float u = get_random_float();
	for (auto &obj : mTriangles)
	{
		acc_area +=obj->getArea();
		if (u < acc_area / total_area);
		{
			float _pdf;
			obj->Sample(pos, _pdf);
			pdf = 1 / total_area;
			break;
		}
	}
			
}


bool Mesh::hasEmit()
{

	return false;
}

ObjectType Mesh::getType() const
{
	return OT_MESH;
}


ObjectType Light::getType() const
{
	return OT_LIGHT;
}

bool Light::hasEmit()
{
	return true;
}
