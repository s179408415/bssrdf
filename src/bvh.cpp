#include "bvh.h"
#include "mesh.h"
#include <algorithm>
//inline bool operator<(const Vector3f& v1, const Vector3f& v2)
//{
//	return 
//}
BVH * BVH::CreateBVH(vector<Object*> objects)
{
	size_t size = objects.size();
	if (size == 0)
		return nullptr;
	AABB aabb;
	for (int i = 0; i < size; i++)
	{
		auto bounds = objects[i]->getBounds();
		aabb += bounds;
		objects[i]->centerid = bounds.GetCenter();
	}
	BVH *bvh = new BVH();
	if (size == 1)
	{
		bvh->mAABB = aabb;
		bvh->mLeft = nullptr;
		bvh->mRight = nullptr;
		bvh->mObject = objects[0];
	}
	else if (size == 2)
	{
		bvh->mLeft = CreateBVH({ objects[0] });
		bvh->mRight = CreateBVH({ objects[1] });
		bvh->mAABB = bvh->mLeft->mAABB;
		bvh->mAABB += bvh->mRight->mAABB;
	}
	else
	{
		Vector3f extent = aabb.GetExtent();
		float Max = max(extent.x(), max(extent.y(), extent.z()));
		int dim = Max == extent.x() ? 0 : (Max == extent.y() ? 1 : 2);
		// everytime sort time killer
		// doesn't needed, so use partial sort, etc std::partition,  std::nth_element,
		// std::partition ,first need known value, the sort is less value on left, bigger value on right, doesn't blance
		// so use std::nth_element make blance

		auto midding = objects.begin() + objects.size() / 2;
		switch (dim)
		{
		case 0:
// 			sort(objects.begin(), objects.end(), [](auto o1, auto o2) {
// 				return o1->getBounds().GetCenter().x() < o2->getBounds().GetCenter().x();
// 				});
			nth_element(objects.begin(), midding, objects.end(), [](auto o1, auto o2) {
				return o1->centerid.x() < o2->centerid.x();
				});
			break;
		case 1:
// 			sort(objects.begin(), objects.end(), [](auto o1, auto o2) {
// 				return o1->getBounds().GetCenter().y() < o2->getBounds().GetCenter().y();
// 				});
			nth_element(objects.begin(), midding, objects.end(), [](auto o1, auto o2) {
				return o1->centerid.y() < o2->centerid.y();
				});
			break;
		case 2:
// 			sort(objects.begin(), objects.end(), [](auto o1, auto o2) {
// 				return o1->getBounds().GetCenter().z() < o2->getBounds().GetCenter().z();
// 				});
			nth_element(objects.begin(), midding, objects.end(), [](auto o1, auto o2) {
				return o1->centerid.z() < o2->centerid.z();
				});
			break;
		default:
			break;
		}

		auto leftObjs = vector<Object*>(objects.begin(), midding);
		auto rightObjs = vector<Object*>(midding, objects.end());
		assert(objects.size() == leftObjs.size() + rightObjs.size());
		bvh->mLeft=CreateBVH(leftObjs);
		bvh->mRight=CreateBVH(rightObjs);
		bvh->mAABB = bvh->mLeft->mAABB;
		bvh->mAABB += bvh->mRight->mAABB;
	}
	return bvh;
}

Intersection BVH::GetIntersection(const Ray& ray) const
{
	Intersection intersect, intersectl, intersectr;
	array<int, 3> dirIsNeg;
	dirIsNeg[0] = int(ray.direction.x() >= 0);
	dirIsNeg[1] = int(ray.direction.y() >= 0);
	dirIsNeg[2] = int(ray.direction.z() >= 0);
	if (!mAABB.IntersectP(ray, ray.direction_inv, dirIsNeg))
		return intersect;
	if (mLeft == nullptr && mRight == nullptr)
	{
		intersect = mObject->getIntersection(ray);
		return intersect;
	}
	intersectl = mLeft->GetIntersection(ray);
	intersectr = mRight->GetIntersection(ray);
	return intersectl.distance < intersectr.distance ? intersectl : intersectr;
}

bool AABB::IntersectP(const Ray& ray, const Vector3f& invDir, const std::array<int, 3>& dirisNeg) const
{
	const auto& origin = ray.origin;
	float ten = -std::numeric_limits<float>::infinity();
	float tex = std::numeric_limits<float>::infinity();
	for (int i = 0; i < 3; i++)
	{
		float min = (Min(i) - origin(i)) * invDir(i);
		float max = (Max(i) - origin(i)) * invDir(i);
		if (!dirisNeg[i])
			std::swap(min, max);
		ten = std::max(min, ten);
		tex = std::min(max, tex);
	}
	return ten <= tex && tex >= 0;
}
