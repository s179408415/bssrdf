#pragma once
#include "assimp/aabb.h"
#include <vector>
#include <Eigen/Eigen>
#include <utility>
using namespace std;
using namespace Eigen;

class Object;

struct Ray {
	//Destination = origin + t*direction
	Vector3f origin;
	Vector3f direction, direction_inv;
	double t;//transportation time,
	double t_min, t_max;
	
	Ray(const Vector3f& ori, const Vector3f& dir, const double _t = 0.0) : origin(ori), direction(dir), t(_t) {
		direction_inv = Vector3f(1. / direction.x(), 1. / direction.y(), 1. / direction.z());
		t_min = 0.0;
		t_max = DBL_MAX; // std::numeric_limits<double>::max();

	}

	Vector3f operator()(double t) const { return origin + direction * t; }

	friend std::ostream& operator<<(std::ostream& os, const Ray& r) {
		os << "[origin:=" << r.origin << ", direction=" << r.direction << ", time=" << r.t << "]\n";
		return os;
	}
};

struct AABB
{
	Vector3f Min;
	Vector3f Max;
	bool bInit = false;
	AABB() {}
	AABB(Vector3f InMin, Vector3f InMax) : Min(InMin), Max(InMax), bInit(true) {}
	AABB(const vector<Vector3f>& vertices) {
		if (vertices.size()>0)
		{
			Min = Vector3f(FLT_MAX, FLT_MAX, FLT_MAX);
			Max = Vector3f(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			for (auto v : vertices)
			{
				for (int i=0;i<3;i++)
				{
					Min(i) = min(Min(i), v(i));
					Max(i) = max(Max(i), v(i));
				}
			}
			bInit = true;
		}
	}

	Vector3f GetExtent() const
	{
		return Max - Min;
	}
	Vector3f GetCenter() const
	{
		return (Min + Max) / 2;
	}
	bool IntersectP(const Ray& ray, const Vector3f& invDir, const std::array<int, 3>& dirisNeg) const;

	AABB &operator+=(const AABB& aabb)
	{
		if (aabb.bInit)
		{
			if (!bInit)
			{
				Min = aabb.Min;
				Max = aabb.Max;
				bInit = true;
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					Min(i) = min(Min(i), aabb.Min(i));
					Max(i) = max(Max(i), aabb.Max(i));
				}

			}
		}
		return *this;
	}
};
struct Intersection
{
	Intersection() {
		happened = false;
		coords = Vector3f();
		normal = Vector3f();
		distance = DBL_MAX;// std::numeric_limits<double>::max();
		obj = nullptr;
		bEmit = false;
		color = Vector3f(0.f, 0.f, 0.f);
		emit = Vector3f(0.f, 0.f, 0.f);
		//m = nullptr;
	}
	bool happened;
	Vector3f coords;
	Vector3f tcoords;
	Vector3f normal;
	Vector3f emit;
	double distance;
	Object* obj;
	Vector3f color;
	bool bEmit;
	//Material* m;
};

class BVH
{
public:
	~BVH() {
		if (mLeft)
			delete mLeft;
		if (mRight)
			delete mRight;
	}

	static BVH* CreateBVH(vector<Object*> objects);

	Intersection GetIntersection(const Ray& ray) const;
	AABB mAABB;
	Object *mObject = nullptr;
	BVH* mLeft = nullptr;
	BVH* mRight = nullptr;
};
