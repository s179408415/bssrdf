#pragma once
#include <vector>
#include <Eigen/Eigen>
#include "bvh.h"
#include "assimp/aabb.h"
#include <random>

using namespace std;
using namespace Eigen;

inline float get_random_float()
{
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_real_distribution<float> dist(0.f, 1.f); // distribution in range [1, 6]

	return dist(rng);
}

enum ObjectType
{
	OT_TRIANGLE,
	OT_MESH,
	OT_LIGHT
};
class Object
{
public:
	Object() {}
	virtual ~Object() {}
	virtual bool intersect(const Ray& ray) = 0;
	virtual bool intersect(const Ray&, float&, uint32_t&) const = 0;
	virtual Intersection getIntersection(Ray _ray) = 0;
	virtual void getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&, Vector2f&) const = 0;
	virtual Vector3f evalDiffuseColor(const Vector2f&) const = 0;
	virtual AABB getBounds() const = 0;
	virtual float getArea() = 0;
	virtual void Sample(Intersection& pos, float& pdf) = 0;
	virtual bool hasEmit() = 0;
	virtual ObjectType getType() const = 0;
	Vector3f centerid;
};

class Triangle : public Object
{
public:
	virtual bool intersect(const Ray& ray);
	virtual bool intersect(const Ray&, float&, uint32_t&) const;
	virtual Intersection getIntersection(Ray _ray);
	virtual void getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&, Vector2f&) const;
	virtual Vector3f evalDiffuseColor(const Vector2f&) const;
	virtual AABB getBounds() const;
	virtual float getArea();
	virtual void Sample(Intersection& pos, float& pdf);
	virtual bool hasEmit();
	virtual ObjectType getType() const override;
	Vector3f v1, v2, v3;
	Vector3f nv1, nv2, nv3;
	Vector3f color = Vector3f(1.f,1.f,1.f);
	bool bEmit = false;
};

class Mesh : public Object
{
public:
	Mesh(const vector<Vector3f>& vertices, const vector<Vector3f>& normals, const vector<int>& indices, Vector3f color = Vector3f(1.f, 1.f, 1.f), bool bEmit=false);
	virtual ~Mesh();
	virtual bool intersect(const Ray& ray);
	virtual bool intersect(const Ray&, float&, uint32_t&) const;
	virtual Intersection getIntersection(Ray _ray);
	virtual void getSurfaceProperties(const Vector3f&, const Vector3f&, const uint32_t&, const Vector2f&, Vector3f&, Vector2f&) const;
	virtual Vector3f evalDiffuseColor(const Vector2f&) const;
	virtual AABB getBounds() const;
	virtual float getArea();
	virtual void Sample(Intersection& pos, float& pdf);
	virtual bool hasEmit();
	virtual ObjectType getType() const override;

	vector<Object*> mTriangles;
	BVH *mBVH = nullptr;
};

class Light : public Mesh
{
public:
	Light(const vector<Vector3f>& vertices, const vector<Vector3f>& normals, const vector<int>& indices, Vector3f color = Vector3f(1.f, 1.f, 1.f))
		:Mesh(vertices, normals, indices, color, true) {}

	virtual ObjectType getType() const override;
	virtual bool hasEmit();
};