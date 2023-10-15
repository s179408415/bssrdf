#include "Scene.h"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include "assimp/aabb.h"
#include "assimp/material.h"
#include "assimp/postprocess.h"
#include <Eigen/Eigen>
using namespace Eigen;
#define PI 3.1415926f

Scene::~Scene()
{
	for (auto &mesh : mMeshes)
	{
		delete mesh;
	}
	if (mBVH)
		delete mBVH;
}

void Scene::LoadScene(const char* filename)
{
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(filename, aiProcess_Triangulate);
	// build model bvh consume much time
	Traversal(scene, scene->mRootNode);
	// add ground plane 
	//CreateGroundPlane();
	// create area light
	//CreateAreaLight();
	CreateBVH();

}

void Scene::CreateBVH()
{
	mBVH = BVH::CreateBVH(mMeshes);
}

void Scene::Remove(int idx)
{
	delete mMeshes[idx];
	mMeshes.erase(mMeshes.begin()+idx);
}

void Scene::CreatePlane(Vector3f pos, Vector2f size, Vector3f Normal, vector<Vector3f> &vertices, vector<Vector3f>& normals, vector<int> &indices)
{
	vertices.resize(0);
	normals.resize(0);
	indices.resize(0);

	vertices.push_back(Vector3f(-size.x(), 0, -size.y()));
	vertices.push_back(Vector3f(-size.x(), 0, size.y()));
	vertices.push_back(Vector3f(size.x(), 0, -size.y()));
	vertices.push_back(Vector3f(size.x(), 0, size.y()));
	normals.push_back(Normal);
	normals.push_back(Normal);
	normals.push_back(Normal);
	normals.push_back(Normal);
	indices.push_back(0);
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(2);
	indices.push_back(1);
	indices.push_back(3);
	Vector3f vx = Normal;
	vx.x() += 0.1f;
	vx.normalize();
	Vector3f vz = vx.cross(Normal).normalized();
	vx = Normal.cross(vz).normalized();
	Matrix3f trans;
	trans << vx, Normal, vz;
	for (auto &vert : vertices)
	{
		vert = trans * vert + pos;
	}
}
void Scene::CreateGroundPlane()
{
	AABB aabb;
	for (int i = 0; i < mMeshes.size(); i++)
	{
		auto _aabb = mMeshes[i]->getBounds();
		aabb += _aabb;
	}
	Vector3f ext = aabb.GetExtent() / 2;
	Vector3f pos = aabb.GetCenter() - Vector3f(0, ext.y(), 0);
	Vector2f SizeXY = Vector2f(ext.x(), ext.z()) * 10.f;
// 	vector<Vector3f> vertices = {
// 		{pos.x() - SizeXY.x(), pos.y(), pos.z() - SizeXY.y()},
// 		{pos.x() - SizeXY.x(), pos.y(), pos.z() + SizeXY.y()},
// 		{pos.x() + SizeXY.x(), pos.y(), pos.z() - SizeXY.y()},
// 		{pos.x() + SizeXY.x(), pos.y(), pos.z() + SizeXY.y()} };
// 	vector<Vector3f> normals = {
// 		{0,1,0},{0,1,0},{0,1,0},{0,1,0}
// 	};
// 	vector<int> indices = {
// 		0,1,2,2,1,3
// 	};
	vector<Vector3f> vertices;
	vector<Vector3f> normals;
	vector<int> indices;
	CreatePlane(pos, SizeXY, Vector3f(0, 1, 0), vertices, normals, indices);
	mMeshes.push_back(new Mesh(vertices, normals, indices));
}

void Scene::CreateAreaLight()
{
	vector<Vector3f> vertices;
	vector<Vector3f> normals;
	vector<int> indices;
	AABB aabb = mMeshes[0]->getBounds();
// 	for (int i = 0; i < mMeshes.size(); i++)
// 	{
// 		auto _aabb = mMeshes[i]->getBounds();
// 		aabb += _aabb;
// 	}
	Vector3f ext = aabb.GetExtent() / 2;
	

	float radius = ext.norm() / 2;
	Vector3f normal(-1, -0.1, 0.3);
	normal.normalize();
	Vector3f pos = aabb.GetCenter() - normal * radius*1.5f + Vector3f(0.f,0.f,-1.f)*radius*0.1f;
	CreatePlane(pos, Vector2f(radius,radius), normal, vertices, normals, indices);
	mMeshes.push_back(new Light(vertices, normals, indices, Vector3f(1.0f,1.f,1.f)));

// 	normal=Vector3f(1, -0.1, 0);
// 	normal.normalize();
// 	pos = aabb.GetCenter() - normal * radius * 1.5f;
// 
// 	CreatePlane(pos, Vector2f(0.05f, 0.05f), normal, vertices, normals, indices);
// 	mMeshes.push_back(new Light(vertices, normals, indices, Vector3f(0.0f, 0.f, 1.f)));

}


Intersection Scene::intersect(const Ray& ray) const
{
	Intersection isect;
	if (!mBVH)
		return isect;
	isect = mBVH->GetIntersection(ray);
	return isect;
}

std::vector<Light*> Scene::GetAllLights() const
{
	vector<Light*> lights;
	for (auto &mesh : mMeshes)
	{
		if (mesh->getType() == OT_LIGHT)
		{
			lights.push_back(static_cast<Light*>(mesh));
		}
	}
	return lights;
}

void Scene::Traversal(const aiScene* scene, aiNode* Node)
{
	for (int i = 0; i < Node->mNumMeshes; i++)
	{
		CreateMesh(scene->mMeshes[Node->mMeshes[i]]);
	}
	for (int i=0;i<Node->mNumChildren;i++)
	{
		Traversal(scene, Node->mChildren[i]);
	}
}

void Scene::CreateMesh(aiMesh* mesh)
{
	vector<Vector3f> vertices;
	vector<Vector3f> normals;
 	float scale = 0.02f;
 	Matrix3f Rotation, Rotation2;
 	float cosval = cos(PI / 2.f);
 	float sinval = sin(PI / 2.f);
 	Rotation << 1.f, 0, 0,
 		0, cosval, -sinval,
 		0, sinval, cosval;
 	Rotation2 << cosval, 0, -sinval,
 		0, 1, 0,
 		sinval, 0, cosval;
	for (int i=0;i<mesh->mNumVertices;i++)
	{
		Vector3f v;
		v << mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z;
 		v *= scale;
		v =  Rotation * Rotation2 * v;
		v += Vector3f(0.2, 0.3, 0.78);
 
		vertices.push_back(v);
		v <<mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z;
		v = Rotation * Rotation2 * v;

		normals.push_back(v);
	}
	vector<int> indices;
	for (int i = 0; i < mesh->mNumFaces; i++)
	{
		if (mesh->mFaces[i].mNumIndices == 3)
		{
			for (int j=0;j<3;j++)
			{
				indices.push_back(mesh->mFaces[i].mIndices[j]);
			}
		}
	}
	mMeshes.push_back(new Mesh(vertices, normals, indices));
}

