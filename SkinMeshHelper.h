//*****************************************************************************************//
//**                                                                                     **//
//**                   　  　　    SkinMeshHelper                                        **//
//**                                                                                     **//
//*****************************************************************************************//

#ifndef Class_SkinMeshHelper_Header
#define Class_SkinMeshHelper_Header

#include "../FbxLoader/FbxLoader.h"
#include "../CoordTf/CoordTf.h"
#include <memory>
#include <stdexcept>

class SkinMeshHelper {

public:
	static const int MAX_BONES = 256;
	static const int FBX_PCS = 32;

	struct Skin_VERTEX {
		CoordTf::VECTOR3 vPos = {};
		CoordTf::VECTOR3 vNorm = {};
		CoordTf::VECTOR3 vTangent = {};
		CoordTf::VECTOR3 vGeoNorm = {};
		CoordTf::VECTOR2 vTex0 = {};
		CoordTf::VECTOR2 vTex1 = {};
		uint32_t bBoneIndex[4] = {};
		float bBoneWeight[4] = {};
	};

	struct Vertex_M {
		CoordTf::VECTOR3 Pos = {};
		CoordTf::VECTOR3 normal = {};
		CoordTf::VECTOR3 tangent = {};
		CoordTf::VECTOR3 geoNormal = {};
		CoordTf::VECTOR2 tex0 = {};
		CoordTf::VECTOR2 tex1 = {};
	};

	struct Skin_VERTEX_Set {
		Skin_VERTEX** pvVB = nullptr;
		Vertex_M** pvVB_M = nullptr;
		uint32_t*** newIndex = nullptr;
		uint32_t** NumNewIndex = nullptr;
	};

	struct BONE {
		CoordTf::MATRIX mBindPose;//初期ポーズ
		CoordTf::MATRIX mNewPose;//現在のポーズ

		BONE()
		{
			memset(this, 0, sizeof(BONE));
		}
	};

	struct SHADER_GLOBAL_BONES {
		CoordTf::MATRIX mBone[MAX_BONES];
		SHADER_GLOBAL_BONES()
		{
			for (int i = 0; i < MAX_BONES; i++)
			{
				CoordTf::MatrixIdentity(&mBone[i]);
			}
		}
	};

	struct meshCenterPos {
		CoordTf::VECTOR3 pos = {};
		uint32_t bBoneIndex = {};
		float bBoneWeight = {};
	};

	void noUseMeshIndex(int meshIndex);
	void ObjCentering(bool f, int ind);
	void ObjCentering(float x, float y, float z, float thetaZ, float thetaY, float thetaX, int ind);
	void ObjOffset(float x, float y, float z, float thetaZ, float thetaY, float thetaX, int ind);
	void SetConnectStep(int ind, float step);

	void GetFbx(char* szFileName);
	void GetFbxSetBinary(char* byteArray, unsigned int size);

	int32_t getMaxEndframe(int fbxIndex, int InternalAnimationIndex);
	int getNumMesh() { return numMesh; }

	void GetFbxSub(char* szFileName, int ind);
	void GetFbxSubSetBinary(char* byteArray, unsigned int size, int ind);
	void GetBuffer_Sub(int ind, int num_end_frame, float* end_frame);
	void GetBuffer_Sub(int ind, float end_frame);
	void CreateFromFBX_SubAnimation(int ind);

	void setDirectTime(float ti);

private:
	class SkinMesh_sub {
	public:
		FbxLoader* fbxL = nullptr;
		std::unique_ptr<float[]> end_frame = nullptr;
		float current_frame = 0.0f;
		bool centering = false;
		bool offset = false;
		float cx = 0.0f;
		float cy = 0.0f;
		float cz = 0.0f;
		float connect_step;
		CoordTf::MATRIX rotZYX = {};
		int InternalLastAnimationIndex = 0;

		SkinMesh_sub();
		~SkinMesh_sub();
		bool Create(char* szFileName);
		bool CreateSetBinary(char* byteArray, unsigned int size);
	};

	//ボーン
	int* numBone = nullptr;
	int maxNumBoneMeshIndex = 0;
	BONE* m_BoneArray = nullptr;
	char* boneName = nullptr;

	std::unique_ptr<meshCenterPos[]> centerPos = nullptr;

	//FBX
	int numMesh = 0;
	SkinMesh_sub* fbx = nullptr;
	Deformer** m_ppSubAnimationBone = nullptr;//その他アニメーションボーンポインタ配列
	CoordTf::MATRIX* m_pLastBoneMatrix = nullptr;
	int AnimLastInd;
	float BoneConnect;
	CoordTf::MATRIX Axis = {};
	std::unique_ptr<bool[]> noUseMesh = nullptr;
	float directTime = 0.0f;

	void DestroyFBX();
	bool InitFBX(char* szFileName, int p);
	bool InitFBXSetBinary(char* byteArray, unsigned int size, int p);
	void CreateRotMatrix(float thetaZ, float thetaY, float thetaX, int ind);
	void ReadSkinInfo(FbxMeshNode* mesh, Skin_VERTEX* tmpVB, meshCenterPos* centerPos);
	void normalRecalculation(bool lclOn, double** nor, FbxMeshNode* mesh);
	void createAxis();
	void AxisSw(bool axisOn);
	void LclTransformation(FbxMeshNode* mesh, CoordTf::VECTOR3* vec);
	void splitIndex(uint32_t numMaterial, FbxMeshNode* mesh, int meshIndex, Skin_VERTEX_Set& vset);

protected:
	int maxNumBone = 0;
	SHADER_GLOBAL_BONES sgb[2] = {};

	SkinMeshHelper();
	~SkinMeshHelper();

	CoordTf::MATRIX GetCurrentPoseMatrix(int index);
	bool SetNewPoseMatrices(float time, int ind, int InternalAnimationIndex);
	void MatrixMap_Bone(SHADER_GLOBAL_BONES* sbB);

	void getBuffer(int num_end_frame, float* end_frame, bool singleMesh, bool deformer);
	Skin_VERTEX_Set setVertex(bool lclOn, bool axisOn, bool VerCentering);

	int getNumBone(int index);
	FbxLoader* getFbxLoader();
	meshCenterPos getMeshCenterPos(int index);
	bool isNoUseMesh(int index);
};

#endif