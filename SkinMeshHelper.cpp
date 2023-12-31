//*****************************************************************************************//
//**                                                                                     **//
//**                   　  　　    SkinMeshHelper                                        **//
//**                                                                                     **//
//*****************************************************************************************//

#define _CRT_SECURE_NO_WARNINGS
#include "SkinMeshHelper.h"

namespace {
	template<typename TYPE>
	void S_DELETE(TYPE p) { if (p) { delete p;    p = nullptr; } }
	template<typename TYPE>
	void A_DELETE(TYPE p) { if (p) { delete[] p;    p = nullptr; } }

	CoordTf::VECTOR3 normalRecalculation_sub(CoordTf::VECTOR3 N[3]) {
		using namespace CoordTf;
		VECTOR3 vecX = {};
		VECTOR3 vecY = {};
		vecX.as(N[0].x - N[1].x, N[0].y - N[1].y, N[0].z - N[1].z);
		vecY.as(N[0].x - N[2].x, N[0].y - N[2].y, N[0].z - N[2].z);
		VECTOR3 vec = {};
		VectorCross(&vec, &vecX, &vecY);
		return vec;
	}

	void createAxisSub(CoordTf::VECTOR3& v3, int axis, int sign) {
		switch (axis) {
		case 0:
			v3.as((float)sign, 0.0f, 0.0f);
			break;
		case 1:
			v3.as(0.0f, (float)sign, 0.0f);
			break;
		case 2:
			v3.as(0.0f, 0.0f, (float)sign);
			break;
		}
	}
}

SkinMeshHelper::SkinMesh_sub::SkinMesh_sub() {
	fbxL = new FbxLoader();
	CoordTf::MatrixIdentity(&rotZYX);
	connect_step = 300.0f;
	InternalLastAnimationIndex = -1;
}

SkinMeshHelper::SkinMesh_sub::~SkinMesh_sub() {
	S_DELETE(fbxL);
}

bool SkinMeshHelper::SkinMesh_sub::Create(char* szFileName) {
	return fbxL->setFbxFile(szFileName);
}

bool SkinMeshHelper::SkinMesh_sub::CreateSetBinary(char* byteArray, unsigned int size) {
	return fbxL->setBinaryInFbxFile(byteArray, size);
}

SkinMeshHelper::SkinMeshHelper() {
	fbx = new SkinMesh_sub[FBX_PCS];
	BoneConnect = -1.0f;
	AnimLastInd = -1;
	directTime = -1.0f;
}

SkinMeshHelper::~SkinMeshHelper() {
	A_DELETE(numBone);
	A_DELETE(boneName);
	A_DELETE(m_ppSubAnimationBone);
	A_DELETE(m_pLastBoneMatrix);
	A_DELETE(m_BoneArray);
	DestroyFBX();
}

bool SkinMeshHelper::InitFBX(char* szFileName, int p) {
	return fbx[p].Create(szFileName);
}

bool SkinMeshHelper::InitFBXSetBinary(char* byteArray, unsigned int size, int p) {
	return fbx[p].CreateSetBinary(byteArray, size);
}

void SkinMeshHelper::DestroyFBX() {
	A_DELETE(fbx);
}

void SkinMeshHelper::ObjCentering(bool f, int ind) {
	fbx[ind].centering = f;
	fbx[ind].offset = false;
	fbx[ind].cx = fbx[ind].cy = fbx[ind].cz = 0.0f;
}

void SkinMeshHelper::CreateRotMatrix(float thetaZ, float thetaY, float thetaX, int ind) {
	using namespace CoordTf;
	MATRIX rotZ, rotY, rotX, rotZY;
	MatrixIdentity(&fbx[ind].rotZYX);
	MatrixRotationZ(&rotZ, thetaZ);
	MatrixRotationY(&rotY, thetaY);
	MatrixRotationX(&rotX, thetaX);
	MatrixMultiply(&rotZY, &rotZ, &rotY);
	MatrixMultiply(&fbx[ind].rotZYX, &rotZY, &rotX);
}

void SkinMeshHelper::ObjCentering(float x, float y, float z, float thetaZ, float thetaY, float thetaX, int ind) {
	fbx[ind].centering = true;
	fbx[ind].offset = false;
	fbx[ind].cx = x;
	fbx[ind].cy = y;
	fbx[ind].cz = z;
	CreateRotMatrix(thetaZ, thetaY, thetaX, ind);
}

void SkinMeshHelper::ObjOffset(float x, float y, float z, float thetaZ, float thetaY, float thetaX, int ind) {
	fbx[ind].centering = true;
	fbx[ind].offset = true;
	fbx[ind].cx = x;
	fbx[ind].cy = y;
	fbx[ind].cz = z;
	CreateRotMatrix(thetaZ, thetaY, thetaX, ind);
}

void SkinMeshHelper::ReadSkinInfo(FbxMeshNode* mesh, Skin_VERTEX* tmpVB, meshCenterPos* cenPos) {

	cenPos->bBoneWeight = 0.0f;
	int addWeightCnt = 0;
	int numbone = mesh->getNumDeformer();
	//各Boneのウエイト,インデックスを調べ頂点配列に加える
	if (numbone <= 0)return;
	for (int i = 0; i < numbone; i++) {
		Deformer* defo = mesh->getDeformer(i);
		int iNumIndex = defo->getIndicesCount();//このボーンに影響を受ける頂点インデックス数
		int* piIndex = defo->getIndices();     //このボーンに影響を受ける頂点のインデックス配列
		double* pdWeight = defo->getWeights();//このボーンに影響を受ける頂点のウエイト配列

		for (int k = 0; k < iNumIndex; k++) {
			int index = piIndex[k];      //影響を受ける頂点
			double weight = pdWeight[k];//ウエイト
			for (int m = 0; m < 4; m++) {
				//各Bone毎に影響を受ける頂点のウエイトを一番大きい数値に更新していく
				if (weight > tmpVB[index].bBoneWeight[m]) {//調べたウエイトの方が大きい
					tmpVB[index].bBoneIndex[m] = i;//Boneインデックス登録
					cenPos->bBoneIndex = i;//とりあえずどれか
					tmpVB[index].bBoneWeight[m] = (float)weight;//ウエイト登録
					cenPos->bBoneWeight += (float)weight;
					addWeightCnt++;
					break;
				}
			}
		}
	}
	cenPos->bBoneWeight /= (float)addWeightCnt;
	//ウエイト正規化
	for (uint32_t i = 0; i < mesh->getNumVertices(); i++) {
		float we = 0;
		for (int j = 0; j < 4; j++) {
			we += tmpVB[i].bBoneWeight[j];
		}
		float we1 = 1.0f / we;
		for (int j = 0; j < 4; j++) {
			tmpVB[i].bBoneWeight[j] *= we1;
		}
	}
}

CoordTf::MATRIX SkinMeshHelper::GetCurrentPoseMatrix(int index) {
	using namespace CoordTf;
	MATRIX inv;
	MatrixIdentity(&inv);
	MatrixInverse(&inv, &m_BoneArray[index].mBindPose);//FBXのバインドポーズは初期姿勢（絶対座標）
	MATRIX fPose;
	MatrixIdentity(&fPose);
	MatrixMultiply(&fPose, &inv, &m_BoneArray[index].mNewPose);//バインドポーズの逆行列とフレーム姿勢行列をかける
	MATRIX ret;
	MatrixIdentity(&ret);
	MatrixMultiply(&ret, &fPose, &Axis);
	return ret;
}

void SkinMeshHelper::MatrixMap_Bone(SHADER_GLOBAL_BONES* sgb) {

	using namespace CoordTf;
	for (int i = 0; i < maxNumBone; i++)
	{
		MATRIX mat = GetCurrentPoseMatrix(i);
		MatrixTranspose(&mat);
		sgb->mBone[i] = mat;
	}
}

bool SkinMeshHelper::SetNewPoseMatrices(float ti, int ind, int InternalAnimationIndex) {
	if (maxNumBone <= 0)return true;
	if (AnimLastInd == -1)AnimLastInd = ind;//最初に描画するアニメーション番号で初期化
	if (fbx[ind].InternalLastAnimationIndex == -1)fbx[ind].InternalLastAnimationIndex = InternalAnimationIndex;

	bool ind_change = false;
	if (AnimLastInd != ind || fbx[ind].InternalLastAnimationIndex != InternalAnimationIndex) {//アニメーションが切り替わった
		ind_change = true;
		AnimLastInd = ind;
		fbx[ind].InternalLastAnimationIndex = InternalAnimationIndex;
		fbx[ind].current_frame = 0.0f;//アニメーションが切り替わってるのでMatrix更新前にフレームを0に初期化
	}
	bool frame_end = false;
	fbx[ind].current_frame += ti;
	if (directTime >= 0.0f) {
		fbx[ind].current_frame = directTime;
	}
	directTime = -1.0f;

	if (fbx[ind].end_frame[InternalAnimationIndex] <= fbx[ind].current_frame) frame_end = true;

	if (frame_end || ind_change) {
		for (int i = 0; i < maxNumBone; i++) {
			for (int x = 0; x < 4; x++) {
				for (int y = 0; y < 4; y++) {
					m_pLastBoneMatrix[i].m[y][x] = m_BoneArray[i].mNewPose.m[y][x];
				}
			}
		}
		BoneConnect = 0.1f;
	}

	if (BoneConnect != -1.0f)fbx[ind].current_frame = 0.0f;

	int frame = (int)fbx[ind].current_frame;
	Deformer Time;
	int64_t time = Time.getTimeFRAMES60(frame);

	bool subanm = true;
	if (ind <= 0 || ind > FBX_PCS - 1)subanm = false;

	Deformer* defo = nullptr;
	FbxLoader* fbL = fbx[0].fbxL;
	FbxMeshNode* mesh = fbL->getFbxMeshNode(maxNumBoneMeshIndex);
	if (!subanm) {
		defo = mesh->getDeformer(0);
	}
	else {
		defo = m_ppSubAnimationBone[(ind - 1) * maxNumBone];
	}
	defo->EvaluateGlobalTransform(time, InternalAnimationIndex);

	using namespace CoordTf;
	MATRIX mat0_mov;
	MatrixIdentity(&mat0_mov);
	if (fbx[ind].offset) {
		MatrixTranslation(&mat0_mov, fbx[ind].cx, fbx[ind].cy, fbx[ind].cz);
	}
	else {
		MatrixTranslation(&mat0_mov, (float)-defo->getEvaluateGlobalTransform(3, 0) + fbx[ind].cx,
			(float)-defo->getEvaluateGlobalTransform(3, 1) + fbx[ind].cy,
			(float)-defo->getEvaluateGlobalTransform(3, 2) + fbx[ind].cz);
	}

	MATRIX pose;
	for (int i = 0; i < maxNumBone; i++) {
		Deformer* de = nullptr;
		if (!subanm) {
			de = mesh->getDeformer(i);
		}
		else {
			de = m_ppSubAnimationBone[(ind - 1) * maxNumBone + i];
		}
		de->EvaluateGlobalTransform(time, InternalAnimationIndex);

		for (int x = 0; x < 4; x++) {
			for (int y = 0; y < 4; y++) {
				if (fbx[ind].centering) {
					pose.m[y][x] = (float)de->getEvaluateGlobalTransform(y, x);
				}
				else {
					m_BoneArray[i].mNewPose.m[y][x] = (float)de->getEvaluateGlobalTransform(y, x);
				}
			}
		}

		if (fbx[ind].centering) {
			MATRIX tmp;
			MatrixMultiply(&tmp, &fbx[ind].rotZYX, &mat0_mov);
			MatrixMultiply(&m_BoneArray[i].mNewPose, &pose, &tmp);
		}
	}

	if (frame_end)fbx[ind].current_frame = 0.0f;

	if (BoneConnect != -1.0f) {
		if (fbx[ind].connect_step <= 0.0f || BoneConnect > 1.0f)BoneConnect = -1.0f;
		else {
			for (int i = 0; i < maxNumBone; i++) {
				StraightLinear(&m_BoneArray[i].mNewPose, &m_pLastBoneMatrix[i], &m_BoneArray[i].mNewPose, BoneConnect += (ti / fbx[ind].connect_step));
			}
		}
	}
	return frame_end;
}

void SkinMeshHelper::normalRecalculation(bool lclOn, double** nor, FbxMeshNode* mesh) {

	*nor = new double[mesh->getNumPolygonVertices() * 3];
	auto index = mesh->getPolygonVertices();//頂点Index取得(頂点xyzに対してのIndex)
	auto ver = mesh->getVertices();//頂点取得

	CoordTf::VECTOR3 tmpv[3] = {};
	uint32_t indexCnt = 0;
	for (uint32_t i = 0; i < mesh->getNumPolygon(); i++) {
		uint32_t pSize = mesh->getPolygonSize(i);//1ポリゴンでの頂点数
		if (pSize >= 3) {
			for (uint32_t i1 = 0; i1 < 3; i1++) {
				uint32_t ind = indexCnt + i1;
				tmpv[i1].as(
					(float)ver[index[ind] * 3],
					(float)ver[index[ind] * 3 + 1],
					(float)ver[index[ind] * 3 + 2]
				);
				if (lclOn) {
					LclTransformation(mesh, &tmpv[i1]);
				}
			}
			//上記3頂点から法線の方向算出
			CoordTf::VECTOR3 N = normalRecalculation_sub(tmpv);
			for (uint32_t ps = 0; ps < pSize; ps++) {
				uint32_t ind = indexCnt + ps;
				(*nor)[ind * 3] = (double)N.x;
				(*nor)[ind * 3 + 1] = (double)N.y;
				(*nor)[ind * 3 + 2] = (double)N.z;
			}
			indexCnt += pSize;
		}
	}
}

void SkinMeshHelper::createAxis() {
	using namespace CoordTf;
	GlobalSettings gSet = fbx[0].fbxL->getGlobalSettings();
	VECTOR3 upVec = {};
	createAxisSub(upVec, gSet.UpAxis, gSet.UpAxisSign);
	VECTOR3 frontVec = {};
	createAxisSub(frontVec, gSet.FrontAxis, gSet.FrontAxisSign);
	VECTOR3 coordVec = {};
	createAxisSub(coordVec, gSet.CoordAxis, gSet.CoordAxisSign);

	MATRIX scale = {};
	MatrixScaling(&scale,
		(float)gSet.UnitScaleFactor,
		(float)gSet.UnitScaleFactor,
		(float)gSet.UnitScaleFactor);

	Axis._11 = coordVec.x; Axis._12 = coordVec.y; Axis._13 = coordVec.z; Axis._14 = 0.0f;
	Axis._21 = upVec.x;    Axis._22 = upVec.y;    Axis._23 = upVec.z;    Axis._24 = 0.0f;
	Axis._31 = frontVec.x; Axis._32 = frontVec.y; Axis._33 = frontVec.z; Axis._34 = 0.0f;
	Axis._41 = 0.0f;       Axis._42 = 0.0f;       Axis._43 = 0.0f;       Axis._44 = 1.0f;

	Axis = Axis * scale;
}

void SkinMeshHelper::LclTransformation(FbxMeshNode* mesh, CoordTf::VECTOR3* vec) {
	using namespace CoordTf;
	MATRIX mov;
	MATRIX rotZ, rotY, rotX, rotZY, rotZYX;
	MATRIX scale;
	MATRIX scro;
	MATRIX world;

	MatrixScaling(&scale,
		(float)mesh->getLcl().Scaling[0],
		(float)mesh->getLcl().Scaling[1],
		(float)mesh->getLcl().Scaling[2]);

	MatrixRotationZ(&rotZ, (float)mesh->getLcl().Rotation[2]);
	MatrixRotationY(&rotY, (float)mesh->getLcl().Rotation[1]);
	MatrixRotationX(&rotX, (float)mesh->getLcl().Rotation[0]);
	MatrixMultiply(&rotZY, &rotZ, &rotY);
	MatrixMultiply(&rotZYX, &rotZY, &rotX);

	MatrixTranslation(&mov,
		(float)mesh->getLcl().Translation[0],
		(float)mesh->getLcl().Translation[1],
		(float)mesh->getLcl().Translation[2]);

	MatrixMultiply(&scro, &scale, &rotZYX);
	MatrixMultiply(&world, &scro, &mov);
	VectorMatrixMultiply(vec, &world);
}

void SkinMeshHelper::SetConnectStep(int ind, float step) {
	fbx[ind].connect_step = step;
}

bool SkinMeshHelper::GetFbx(char* szFileName) {
	//FBXローダーを初期化
	return InitFBX(szFileName, 0);
}

bool SkinMeshHelper::GetFbxSetBinary(char* byteArray, unsigned int size) {
	//FBXローダーを初期化
	return InitFBXSetBinary(byteArray, size, 0);
}

int32_t SkinMeshHelper::getMaxEndframe(int fbxIndex, int InternalAnimationIndex) {
	FbxLoader* fl = fbx[fbxIndex].fbxL;
	FbxMeshNode* mesh = fl->getFbxMeshNode(0);
	Deformer* defo = nullptr;
	if (fbxIndex == 0) {
		defo = mesh->getDeformer(0);
	}
	else {
		defo = fl->getNoneMeshDeformer(0);
	}
	return defo->getMaxFRAMES60(InternalAnimationIndex);
}

void SkinMeshHelper::AxisSw(bool axisOn) {
	if (axisOn) {
		createAxis();
	}
	else {
		CoordTf::MatrixIdentity(&Axis);
	}
}

Skin_VERTEX_Set SkinMeshHelper::setVertex(bool lclOn, bool axisOn, bool VerCentering) {

	AxisSw(axisOn);

	Skin_VERTEX_Set vset = {};

	vset.pvVB = new Skin_VERTEX * [numMesh];
	vset.pvVB_M = new Vertex_M * [numMesh];
	vset.newIndex = new uint32_t * *[numMesh];
	vset.NumNewIndex = new uint32_t * [numMesh];

	for (int m = 0; m < numMesh; m++) {

		FbxLoader* fbL = fbx[0].fbxL;
		FbxMeshNode* mesh = fbL->getFbxMeshNode(m);//メッシュ毎に処理する

		//分割後のIndex生成
		splitIndex(mesh->getNumMaterial(), mesh, m, vset);

		Skin_VERTEX* tmpVB = new Skin_VERTEX[mesh->getNumVertices()];
		//ボーンウエイト
		ReadSkinInfo(mesh, tmpVB, &centerPos[m]);

		auto index = mesh->getPolygonVertices();//頂点Index取得(頂点xyzに対してのIndex)
		auto ver = mesh->getVertices();//頂点取得
		auto nor = mesh->getAlignedNormal(0);//法線取得

		bool norCreate = false;
		if (!nor) {
			normalRecalculation(lclOn, &nor, mesh);
			norCreate = true;
		}

		double* uv0 = mesh->getAlignedUV(0);//テクスチャUV0
		double* uv1 = nullptr;              //テクスチャUV1

		if (mesh->getNumUVObj() > 1) {
			uv1 = mesh->getAlignedUV(1);
		}
		else {
			uv1 = mesh->getAlignedUV(0);
		}

		Skin_VERTEX* vb = vset.pvVB[m] = new Skin_VERTEX[mesh->getNumPolygonVertices()];
		Vertex_M* vbm = vset.pvVB_M[m] = new Vertex_M[mesh->getNumPolygonVertices()];

		meshCenterPos& cp = centerPos[m];
		CoordTf::VECTOR3& cpp = cp.pos;
		int cppAddCnt = 0;

		for (uint32_t i = 0; i < mesh->getNumPolygonVertices(); i++) {
			//index順で頂点を整列しながら頂点格納
			Skin_VERTEX* v = &vb[i];
			Vertex_M* vm = &vbm[i];
			v->vPos.x = (float)ver[index[i] * 3];
			v->vPos.y = (float)ver[index[i] * 3 + 1];
			v->vPos.z = (float)ver[index[i] * 3 + 2];
			if (lclOn) {
				LclTransformation(mesh, &v->vPos);
			}
			cpp.x += vm->Pos.x = v->vPos.x;
			cpp.y += vm->Pos.y = v->vPos.y;
			cpp.z += vm->Pos.z = v->vPos.z;
			cppAddCnt++;

			int norInd = i;
			int uvInd = i;
			if (mesh->getNormalMappingInformationType(0) &&
				strcmp(mesh->getNormalMappingInformationType(0), "ByPolygonVertex")) {
				norInd = index[i];
			}

			if (mesh->getUVMappingInformationType(0) &&
				strcmp(mesh->getUVMappingInformationType(0), "ByPolygonVertex")) {
				uvInd = index[i];
			}

			vm->normal.x = v->vNorm.x = (float)nor[norInd * 3];
			vm->normal.y = v->vNorm.y = (float)nor[norInd * 3 + 1];
			vm->normal.z = v->vNorm.z = (float)nor[norInd * 3 + 2];
			vm->geoNormal.x = v->vGeoNorm.x = (float)nor[norInd * 3];
			vm->geoNormal.y = v->vGeoNorm.y = (float)nor[norInd * 3 + 1];
			vm->geoNormal.z = v->vGeoNorm.z = (float)nor[norInd * 3 + 2];
			vm->tex0.x = v->vTex0.x = (float)uv0[uvInd * 2];
			vm->tex0.y = v->vTex0.y = 1.0f - (float)uv0[uvInd * 2 + 1];//(1.0f-UV)
			vm->tex1.x = v->vTex1.x = (float)uv1[uvInd * 2];
			vm->tex1.y = v->vTex1.y = 1.0f - (float)uv1[uvInd * 2 + 1];//(1.0f-UV)

			if (numBone[m] > 0) {
				for (int bi = 0; bi < 4; bi++) {
					//ReadSkinInfo(tmpVB)で読み込んだ各パラメータコピー
					v->bBoneIndex[bi] = tmpVB[index[i]].bBoneIndex[bi];
					v->bBoneWeight[bi] = tmpVB[index[i]].bBoneWeight[bi];
				}
			}
		}
		A_DELETE(tmpVB);
		if (norCreate)A_DELETE(nor);

		cpp.x /= (float)cppAddCnt;
		cpp.y /= (float)cppAddCnt;
		cpp.z /= (float)cppAddCnt;
		if (VerCentering) {
			for (uint32_t i = 0; i < mesh->getNumPolygonVertices(); i++) {
				vb[i].vPos.x -= cpp.x;
				vb[i].vPos.y -= cpp.y;
				vb[i].vPos.z -= cpp.z;
				vbm[i].Pos.x -= cpp.x;
				vbm[i].Pos.y -= cpp.y;
				vbm[i].Pos.z -= cpp.z;
			}
			cpp.x = 0.0f;
			cpp.y = 0.0f;
			cpp.z = 0.0f;
		}
	}

	return vset;
}

void SkinMeshHelper::splitIndex(uint32_t numMaterial, FbxMeshNode* mesh, int m, Skin_VERTEX_Set& vset) {

	//ポリゴン分割後のIndex数カウント
	vset.NumNewIndex[m] = new uint32_t[numMaterial];
	uint32_t* numNewIndex = vset.NumNewIndex[m];

	memset(numNewIndex, 0, sizeof(uint32_t) * numMaterial);

	for (uint32_t i = 0; i < mesh->getNumPolygon(); i++) {
		uint32_t currentMatNo = mesh->getMaterialNoOfPolygon(i);
		uint32_t pSize = mesh->getPolygonSize(i);
		if (pSize >= 3) {
			numNewIndex[currentMatNo] += (pSize - 2) * 3;
		}
	}

	//分割後のIndex生成
	vset.newIndex[m] = new uint32_t * [numMaterial];
	uint32_t** newIndex = vset.newIndex[m];
	for (uint32_t ind1 = 0; ind1 < numMaterial; ind1++) {
		if (numNewIndex[ind1] <= 0) {
			newIndex[ind1] = nullptr;
			continue;
		}
		newIndex[ind1] = new uint32_t[numNewIndex[ind1]];
	}
	std::unique_ptr<uint32_t[]> indexCnt;
	indexCnt = std::make_unique<uint32_t[]>(numMaterial);

	int Icnt = 0;
	for (uint32_t i = 0; i < mesh->getNumPolygon(); i++) {
		uint32_t currentMatNo = mesh->getMaterialNoOfPolygon(i);
		uint32_t pSize = mesh->getPolygonSize(i);
		if (pSize >= 3) {
			for (uint32_t ps = 0; ps < pSize - 2; ps++) {
				newIndex[currentMatNo][indexCnt[currentMatNo]++] = Icnt;
				newIndex[currentMatNo][indexCnt[currentMatNo]++] = Icnt + 1 + ps;
				newIndex[currentMatNo][indexCnt[currentMatNo]++] = Icnt + 2 + ps;
			}
			Icnt += pSize;
		}
	}
}

void SkinMeshHelper::getBuffer(int num_end_frame, float* end_frame, bool singleMesh, bool deformer) {

	using namespace CoordTf;
	fbx[0].end_frame = std::make_unique<float[]>(num_end_frame);
	memcpy(fbx[0].end_frame.get(), end_frame, num_end_frame * sizeof(float));
	FbxLoader* fbL = fbx[0].fbxL;
	if (singleMesh)fbL->createFbxSingleMeshNode();
	numMesh = fbL->getNumFbxMeshNode();
	noUseMesh = std::make_unique<bool[]>(numMesh);
	centerPos = std::make_unique<meshCenterPos[]>(numMesh);
	numBone = new int[numMesh];
	for (int i = 0; i < numMesh; i++) {
		noUseMesh[i] = false;
		if (deformer) {
			numBone[i] = fbL->getFbxMeshNode(i)->getNumDeformer();
			if (maxNumBone < numBone[i]) {
				maxNumBone = numBone[i];
				maxNumBoneMeshIndex = i;
			}
		}
		else
			numBone[i] = 0;
	}

	if (maxNumBone > 0) {
		boneName = new char[maxNumBone * 255];
		m_BoneArray = new BONE[maxNumBone];
		m_pLastBoneMatrix = new MATRIX[maxNumBone];

		FbxMeshNode* mesh = fbL->getFbxMeshNode(maxNumBoneMeshIndex);
		for (int i = 0; i < maxNumBone; i++) {
			Deformer* defo = mesh->getDeformer(i);
			const char* name = defo->getName();
			strcpy(&boneName[i * 255], name);//ボーンの名前保持

			//初期姿勢行列読み込み
			//GetCurrentPoseMatrixで使う
			for (int x = 0; x < 4; x++) {
				for (int y = 0; y < 4; y++) {
					m_BoneArray[i].mBindPose.m[y][x] = (float)defo->getTransformLinkMatrix(y, x);
				}
			}
		}
	}
}

void SkinMeshHelper::noUseMeshIndex(int meshIndex) {
	noUseMesh[meshIndex] = true;
}

bool SkinMeshHelper::GetFbxSub(char* szFileName, int ind) {
	if (ind <= 0)return false;
	return InitFBX(szFileName, ind);
}

bool SkinMeshHelper::GetFbxSubSetBinary(char* byteArray, unsigned int size, int ind) {
	if (ind <= 0)return false;
	return InitFBXSetBinary(byteArray, size, ind);
}

bool SkinMeshHelper::GetBuffer_Sub(int ind, float end_frame) {
	float ef[1] = { end_frame };
	return GetBuffer_Sub(ind, 1, ef);
}

bool SkinMeshHelper::GetBuffer_Sub(int ind, int num_end_frame, float* end_frame) {

	fbx[ind].end_frame = std::make_unique<float[]>(num_end_frame);
	memcpy(fbx[ind].end_frame.get(), end_frame, num_end_frame * sizeof(float));

	int BoneNum = fbx[ind].fbxL->getNumNoneMeshDeformer();
	if (BoneNum == 0) return false;

	if (!m_ppSubAnimationBone) {
		m_ppSubAnimationBone = new Deformer * [(FBX_PCS - 1) * maxNumBone];
	}
	return true;
}

void SkinMeshHelper::CreateFromFBX_SubAnimation(int ind) {
	int loopind = 0;
	int searchCount = 0;
	int name2Num = 0;
	while (loopind < maxNumBone) {
		int sa_ind = (ind - 1) * maxNumBone + loopind;
		m_ppSubAnimationBone[sa_ind] = fbx[ind].fbxL->getNoneMeshDeformer(searchCount);
		searchCount++;
		const char* name = m_ppSubAnimationBone[sa_ind]->getName();
		if (!strncmp("Skeleton", name, 8))continue;
		char* name2 = &boneName[loopind * 255];//各Bone名の先頭アドレスを渡す
		//Bone名に空白が含まれている場合最後の空白以降の文字から終端までの文字を比較する為,
		//終端文字までポインタを進め, 終端から検査して空白位置手前まで比較する
		while (*name != '\0')name++;//終端文字までポインタを進める
		//終端文字までポインタを進める, 空白が含まれない文字の場合もあるので文字数カウントし,
		//文字数で比較完了を判断する
		while (*name2 != '\0') { name2++; name2Num++; }
		while (*(--name) == *(--name2) && *name2 != ' ' && (--name2Num) > 0);
		if (*name2 != ' ' && name2Num > 0) { name2Num = 0; continue; }
		name2Num = 0;
		loopind++;
	}
}

void SkinMeshHelper::setDirectTime(float ti) {
	directTime = ti;
}