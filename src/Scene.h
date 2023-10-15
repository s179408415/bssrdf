#pragma once



#include "Scene.h"
#include "bvh.h"
#include "assimp/scene.h"
#include "mesh.h"
#include <vector>
using namespace std;

class Scene
{
public:
	~Scene();
	void LoadScene(const char* filename);

	void CreateBVH();

	void Remove(int idx);

	void CreatePlane(Vector3f pos, Vector2f size, Vector3f Normal, vector<Vector3f>& vertices, vector<Vector3f>& normals, vector<int>& indices);
	void CreateGroundPlane();

	Intersection intersect(const Ray& ray) const;
	Object* GetObject(int i) const { return mMeshes[i]; }
	vector<Light*> GetAllLights() const;
private:
	BVH *mBVH = nullptr;
	vector<Object*> mMeshes;
	void Traversal(const aiScene* scene, aiNode* Node);
	void CreateMesh(aiMesh* mesh);
	void CreateAreaLight();
};

