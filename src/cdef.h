/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/
typedef long dint64;
typedef unsigned long duint64;
typedef int dint32;
typedef unsigned int duint32;
typedef short dint16;
typedef unsigned short duint16;
typedef signed char dint8;
typedef unsigned char duint8;
typedef dint64 dintptr;
typedef duint64 duintptr;
typedef dint64 ddiffint;
typedef duint64 dsizeint;
typedef void dMessageFunction (int errnum, const char *msg, va_list ap);
void dSetErrorHandler (dMessageFunction *fn);
void dSetDebugHandler (dMessageFunction *fn);
void dSetMessageHandler (dMessageFunction *fn);
dMessageFunction *dGetErrorHandler(void);
dMessageFunction *dGetDebugHandler(void);
dMessageFunction *dGetMessageHandler(void);
void dError (int num, const char *msg, ...);
void dDebug (int num, const char *msg, ...);
void dMessage (int num, const char *msg, ...);
typedef double dReal;
typedef duint32 dTriIndex;
typedef enum {
    dSA__MIN, dSA_X = dSA__MIN, dSA_Y, dSA_Z, dSA__MAX,
} dSpaceAxis;
typedef enum {
    dMD__MIN, dMD_LINEAR = dMD__MIN, dMD_ANGULAR, dMD__MAX,
} dMotionDynamics;
typedef enum {
    dDA__MIN, dDA__L_MIN = dDA__MIN + dMD_LINEAR * dSA__MAX, dDA_LX = dDA__L_MIN + dSA_X, dDA_LY = dDA__L_MIN + dSA_Y, dDA_LZ = dDA__L_MIN + dSA_Z, dDA__L_MAX = dDA__L_MIN + dSA__MAX, dDA__A_MIN = dDA__MIN + dMD_ANGULAR * dSA__MAX, dDA_AX = dDA__A_MIN + dSA_X, dDA_AY = dDA__A_MIN + dSA_Y, dDA_AZ = dDA__A_MIN + dSA_Z, dDA__A_MAX = dDA__A_MIN + dSA__MAX, dDA__MAX = dDA__MIN + dMD__MAX * dSA__MAX,
} dDynamicsAxis;
typedef enum {
    dV3E__MIN, dV3E__AXES_MIN = dV3E__MIN, dV3E_X = dV3E__AXES_MIN + dSA_X, dV3E_Y = dV3E__AXES_MIN + dSA_Y, dV3E_Z = dV3E__AXES_MIN + dSA_Z, dV3E__AXES_MAX = dV3E__AXES_MIN + dSA__MAX, dV3E_PAD = dV3E__AXES_MAX, dV3E__MAX, dV3E__AXES_COUNT = dV3E__AXES_MAX - dV3E__AXES_MIN,
} dVec3Element;
typedef enum {
    dV4E__MIN, dV4E_X = dV4E__MIN + dSA_X, dV4E_Y = dV4E__MIN + dSA_Y, dV4E_Z = dV4E__MIN + dSA_Z, dV4E_O = dV4E__MIN + dSA__MAX, dV4E__MAX,
} dVec4Element;
typedef enum {
    dM3E__MIN, dM3E__X_MIN = dM3E__MIN + dSA_X * dV3E__MAX, dM3E__X_AXES_MIN = dM3E__X_MIN + dV3E__AXES_MIN, dM3E_XX = dM3E__X_MIN + dV3E_X, dM3E_XY = dM3E__X_MIN + dV3E_Y, dM3E_XZ = dM3E__X_MIN + dV3E_Z, dM3E__X_AXES_MAX = dM3E__X_MIN + dV3E__AXES_MAX, dM3E_XPAD = dM3E__X_MIN + dV3E_PAD, dM3E__X_MAX = dM3E__X_MIN + dV3E__MAX, dM3E__Y_MIN = dM3E__MIN + dSA_Y * dV3E__MAX, dM3E__Y_AXES_MIN = dM3E__Y_MIN + dV3E__AXES_MIN, dM3E_YX = dM3E__Y_MIN + dV3E_X, dM3E_YY = dM3E__Y_MIN + dV3E_Y, dM3E_YZ = dM3E__Y_MIN + dV3E_Z, dM3E__Y_AXES_MAX = dM3E__Y_MIN + dV3E__AXES_MAX, dM3E_YPAD = dM3E__Y_MIN + dV3E_PAD, dM3E__Y_MAX = dM3E__Y_MIN + dV3E__MAX, dM3E__Z_MIN = dM3E__MIN + dSA_Z * dV3E__MAX, dM3E__Z_AXES_MIN = dM3E__Z_MIN + dV3E__AXES_MIN, dM3E_ZX = dM3E__Z_MIN + dV3E_X, dM3E_ZY = dM3E__Z_MIN + dV3E_Y, dM3E_ZZ = dM3E__Z_MIN + dV3E_Z, dM3E__Z_AXES_MAX = dM3E__Z_MIN + dV3E__AXES_MAX, dM3E_ZPAD = dM3E__Z_MIN + dV3E_PAD, dM3E__Z_MAX = dM3E__Z_MIN + dV3E__MAX, dM3E__MAX = dM3E__MIN + dSA__MAX * dV3E__MAX,
} dMat3Element;
typedef enum {
    dM4E__MIN, dM4E__X_MIN = dM4E__MIN + dV4E_X * dV4E__MAX, dM4E_XX = dM4E__X_MIN + dV4E_X, dM4E_XY = dM4E__X_MIN + dV4E_Y, dM4E_XZ = dM4E__X_MIN + dV4E_Z, dM4E_XO = dM4E__X_MIN + dV4E_O, dM4E__X_MAX = dM4E__X_MIN + dV4E__MAX, dM4E__Y_MIN = dM4E__MIN + dV4E_Y * dV4E__MAX, dM4E_YX = dM4E__Y_MIN + dV4E_X, dM4E_YY = dM4E__Y_MIN + dV4E_Y, dM4E_YZ = dM4E__Y_MIN + dV4E_Z, dM4E_YO = dM4E__Y_MIN + dV4E_O, dM4E__Y_MAX = dM4E__Y_MIN + dV4E__MAX, dM4E__Z_MIN = dM4E__MIN + dV4E_Z * dV4E__MAX, dM4E_ZX = dM4E__Z_MIN + dV4E_X, dM4E_ZY = dM4E__Z_MIN + dV4E_Y, dM4E_ZZ = dM4E__Z_MIN + dV4E_Z, dM4E_ZO = dM4E__Z_MIN + dV4E_O, dM4E__Z_MAX = dM4E__Z_MIN + dV4E__MAX, dM4E__O_MIN = dM4E__MIN + dV4E_O * dV4E__MAX, dM4E_OX = dM4E__O_MIN + dV4E_X, dM4E_OY = dM4E__O_MIN + dV4E_Y, dM4E_OZ = dM4E__O_MIN + dV4E_Z, dM4E_OO = dM4E__O_MIN + dV4E_O, dM4E__O_MAX = dM4E__O_MIN + dV4E__MAX, dM4E__MAX = dM4E__MIN + dV4E__MAX * dV4E__MAX,
} dMat4Element;
typedef enum {
    dQUE__MIN, dQUE_R = dQUE__MIN, dQUE__AXIS_MIN, dQUE_I = dQUE__AXIS_MIN + dSA_X, dQUE_J = dQUE__AXIS_MIN + dSA_Y, dQUE_K = dQUE__AXIS_MIN + dSA_Z, dQUE__AXIS_MAX = dQUE__AXIS_MIN + dSA__MAX, dQUE__MAX = dQUE__AXIS_MAX,
} dQuatElement;
typedef dReal dVector3[dV3E__MAX];
typedef dReal dVector4[dV4E__MAX];
typedef dReal dMatrix3[dM3E__MAX];
typedef dReal dMatrix4[dM4E__MAX];
typedef dReal dMatrix6[(dMD__MAX * dV3E__MAX) * (dMD__MAX * dSA__MAX)];
typedef dReal dQuaternion[dQUE__MAX];
struct dxWorld; 
struct dxSpace; 
struct dxBody; 
struct dxGeom; 
struct dxJoint; 
struct dxJointGroup;
typedef struct dxWorld *dWorldID;
typedef struct dxSpace *dSpaceID;
typedef struct dxBody *dBodyID;
typedef struct dxGeom *dGeomID;
typedef struct dxJoint *dJointID;
typedef struct dxJointGroup *dJointGroupID;
enum {
  d_ERR_UNKNOWN = 0, 
  d_ERR_IASSERT, 
  d_ERR_UASSERT, 
  d_ERR_LCP 
};
typedef enum {
  dJointTypeNone = 0, 
  dJointTypeBall, dJointTypeHinge, dJointTypeSlider, dJointTypeContact, dJointTypeUniversal, dJointTypeHinge2, dJointTypeFixed, dJointTypeNull, dJointTypeAMotor, dJointTypeLMotor, dJointTypePlane2D, dJointTypePR, dJointTypePU, dJointTypePiston, dJointTypeDBall, dJointTypeDHinge, dJointTypeTransmission,
} dJointType;
enum {
  dParamLoStop = 0, dParamHiStop, dParamVel, dParamLoVel, dParamHiVel, dParamFMax, dParamFudgeFactor, dParamBounce, dParamCFM, dParamStopERP, dParamStopCFM, dParamSuspensionERP, dParamSuspensionCFM, dParamERP, dParamsInGroup, 
  dParamGroup1 = 0x000, dParamLoStop1 = 0x000, dParamHiStop1, dParamVel1, dParamLoVel1, dParamHiVel1, dParamFMax1, dParamFudgeFactor1, dParamBounce1, dParamCFM1, dParamStopERP1, dParamStopCFM1, dParamSuspensionERP1, dParamSuspensionCFM1, dParamERP1, dParamGroup2 = 0x100, dParamLoStop2 = 0x100, dParamHiStop2, dParamVel2, dParamLoVel2, dParamHiVel2, dParamFMax2, dParamFudgeFactor2, dParamBounce2, dParamCFM2, dParamStopERP2, dParamStopCFM2, dParamSuspensionERP2, dParamSuspensionCFM2, dParamERP2, dParamGroup3 = 0x200, dParamLoStop3 = 0x200, dParamHiStop3, dParamVel3, dParamLoVel3, dParamHiVel3, dParamFMax3, dParamFudgeFactor3, dParamBounce3, dParamCFM3, dParamStopERP3, dParamStopCFM3, dParamSuspensionERP3, dParamSuspensionCFM3, dParamERP3, 
  dParamGroup=0x100
};
enum {
  dAMotorUser = 0, dAMotorEuler = 1
};
enum {
  dTransmissionParallelAxes = 0, dTransmissionIntersectingAxes = 1, dTransmissionChainDrive = 2
};
typedef struct dJointFeedback {
  dVector3 f1; 
  dVector3 t1; 
  dVector3 f2; 
  dVector3 t2; 
} dJointFeedback;
void dGeomMoved (dGeomID);
dGeomID dGeomGetBodyNext (dGeomID);
const char* dGetConfiguration (void);
int dCheckConfiguration( const char* token );
enum dInitODEFlags {
    dInitFlagManualThreadCleanup = 0x00000001 
};
void dInitODE(void);
int dInitODE2(unsigned int uiInitFlags);
enum dAllocateODEDataFlags {
    dAllocateFlagBasicData = 0, 
    dAllocateFlagCollisionData = 0x00000001, 
    dAllocateMaskAll = ~0 
};
int dAllocateODEDataForThread(unsigned int uiAllocateFlags);
void dCleanupODEAllDataForThread();
void dCloseODE(void);
enum {
  dContactMu2 = 0x001, 
  dContactAxisDep = 0x001, 
  dContactFDir1 = 0x002, 
  dContactBounce = 0x004, 
  dContactSoftERP = 0x008, 
  dContactSoftCFM = 0x010, 
  dContactMotion1 = 0x020, 
  dContactMotion2 = 0x040, dContactMotionN = 0x080, dContactSlip1 = 0x100, 
  dContactSlip2 = 0x200, dContactRolling = 0x400, 
  dContactApprox0 = 0x0000, dContactApprox1_1 = 0x1000, dContactApprox1_2 = 0x2000, dContactApprox1_N = 0x4000, 
  dContactApprox1 = 0x7000
};
typedef struct dSurfaceParameters {
  int mode;
  dReal mu;
  dReal mu2;
  dReal rho; 
  dReal rho2;
  dReal rhoN; 
  dReal bounce; 
  dReal bounce_vel; 
  dReal soft_erp;
  dReal soft_cfm;
  dReal motion1,motion2,motionN;
  dReal slip1,slip2;
} dSurfaceParameters;
typedef struct dContactGeom {
    dVector3 pos; 
    dVector3 normal; 
    dReal depth; 
    dGeomID g1,g2; 
    int side1,side2; 
} dContactGeom;
typedef struct dContact {
  dSurfaceParameters surface;
  dContactGeom geom;
  dVector3 fdir1;
} dContact;
typedef void * dAllocFunction (dsizeint size);
typedef void * dReallocFunction (void *ptr, dsizeint oldsize, dsizeint newsize);
typedef void dFreeFunction (void *ptr, dsizeint size);
void dSetAllocHandler (dAllocFunction *fn);
void dSetReallocHandler (dReallocFunction *fn);
void dSetFreeHandler (dFreeFunction *fn);
dAllocFunction *dGetAllocHandler (void);
dReallocFunction *dGetReallocHandler (void);
dFreeFunction *dGetFreeHandler (void);
void * dAlloc (dsizeint size);
void * dRealloc (void *ptr, dsizeint oldsize, dsizeint newsize);
void dFree (void *ptr, dsizeint size);
int dSafeNormalize3 (dVector3 a);
int dSafeNormalize4 (dVector4 a);
void dNormalize3 (dVector3 a);
void dNormalize4 (dVector4 a);
void dPlaneSpace (const dVector3 n, dVector3 p, dVector3 q);
int dOrthogonalizeR(dMatrix3 m);
void dSetZero (dReal *a, int n);
void dSetValue (dReal *a, int n, dReal value);
dReal dDot (const dReal *a, const dReal *b, int n);
void dMultiply0 (dReal *A, const dReal *B, const dReal *C, int p,int q,int r);
void dMultiply1 (dReal *A, const dReal *B, const dReal *C, int p,int q,int r);
void dMultiply2 (dReal *A, const dReal *B, const dReal *C, int p,int q,int r);
int dFactorCholesky (dReal *A, int n);
void dSolveCholesky (const dReal *L, dReal *b, int n);
int dInvertPDMatrix (const dReal *A, dReal *Ainv, int n);
int dIsPositiveDefinite (const dReal *A, int n);
void dFactorLDLT (dReal *A, dReal *d, int n, int nskip);
void dSolveL1 (const dReal *L, dReal *b, int n, int nskip);
void dSolveL1T (const dReal *L, dReal *b, int n, int nskip);
void dScaleVector (dReal *a, const dReal *d, int n);
void dVectorScale (dReal *a, const dReal *d, int n);
void dSolveLDLT (const dReal *L, const dReal *d, dReal *b, int n, int nskip);
void dLDLTAddTL (dReal *L, dReal *d, const dReal *a, int n, int nskip);
void dLDLTRemove (dReal **A, const int *p, dReal *L, dReal *d, int n1, int n2, int r, int nskip);
void dRemoveRowCol (dReal *A, int n, int nskip, int r);
struct dxThreadingImplementation;
typedef struct dxThreadingImplementation *dThreadingImplementationID;
typedef unsigned dmutexindex_t;
struct dxMutexGroup;
typedef struct dxMutexGroup *dMutexGroupID;
typedef dMutexGroupID dMutexGroupAllocFunction (dThreadingImplementationID impl, dmutexindex_t Mutex_count, const char *const *Mutex_names_ptr);
typedef void dMutexGroupFreeFunction (dThreadingImplementationID impl, dMutexGroupID mutex_group);
typedef void dMutexGroupMutexLockFunction (dThreadingImplementationID impl, dMutexGroupID mutex_group, dmutexindex_t mutex_index);
typedef void dMutexGroupMutexUnlockFunction (dThreadingImplementationID impl, dMutexGroupID mutex_group, dmutexindex_t mutex_index);
struct dxCallReleasee;
typedef struct dxCallReleasee *dCallReleaseeID;
struct dxCallWait;
typedef struct dxCallWait *dCallWaitID;
typedef dsizeint ddependencycount_t;
typedef ddiffint ddependencychange_t;
typedef dsizeint dcallindex_t;
typedef int dThreadedCallFunction(void *call_context, dcallindex_t instance_index, dCallReleaseeID this_releasee);
typedef struct dxThreadedWaitTime
{
  int64_t wait_sec;
  unsigned long wait_nsec;
} dThreadedWaitTime;
typedef dCallWaitID dThreadedCallWaitAllocFunction(dThreadingImplementationID impl);
typedef void dThreadedCallWaitResetFunction(dThreadingImplementationID impl, dCallWaitID call_wait);
typedef void dThreadedCallWaitFreeFunction(dThreadingImplementationID impl, dCallWaitID call_wait);
typedef void dThreadedCallPostFunction(dThreadingImplementationID impl, int *out_summary_fault, dCallReleaseeID *out_post_releasee, ddependencycount_t dependencies_count, dCallReleaseeID dependent_releasee, dCallWaitID call_wait, dThreadedCallFunction *call_func, void *call_context, dcallindex_t instance_index, const char *call_name);
typedef void dThreadedCallDependenciesCountAlterFunction(dThreadingImplementationID impl, dCallReleaseeID target_releasee, ddependencychange_t dependencies_count_change);
typedef void dThreadedCallWaitFunction(dThreadingImplementationID impl, int *out_wait_status, dCallWaitID call_wait, const dThreadedWaitTime *timeout_time_ptr, const char *wait_name);
typedef unsigned dThreadingImplThreadCountRetrieveFunction(dThreadingImplementationID impl);
typedef int dThreadingImplResourcesForCallsPreallocateFunction(dThreadingImplementationID impl, ddependencycount_t max_simultaneous_calls_estimate);
typedef struct dxThreadingFunctionsInfo
{
  unsigned struct_size;
  dMutexGroupAllocFunction *alloc_mutex_group;
  dMutexGroupFreeFunction *free_mutex_group;
  dMutexGroupMutexLockFunction *lock_group_mutex;
  dMutexGroupMutexUnlockFunction *unlock_group_mutex;
  dThreadedCallWaitAllocFunction *alloc_call_wait;
  dThreadedCallWaitResetFunction *reset_call_wait;
  dThreadedCallWaitFreeFunction *free_call_wait;
  dThreadedCallPostFunction *post_call;
  dThreadedCallDependenciesCountAlterFunction *alter_call_dependencies_count;
  dThreadedCallWaitFunction *wait_call;
  dThreadingImplThreadCountRetrieveFunction *retrieve_thread_count;
  dThreadingImplResourcesForCallsPreallocateFunction *preallocate_resources_for_calls;
} dThreadingFunctionsInfo;
struct dxCooperative;
struct dxResourceRequirements;
struct dxResourceContainer;
typedef struct dxCooperative *dCooperativeID;
typedef struct dxResourceRequirements *dResourceRequirementsID;
typedef struct dxResourceContainer *dResourceContainerID;
dCooperativeID dCooperativeCreate(const dThreadingFunctionsInfo *functionInfo, dThreadingImplementationID threadingImpl);
void dCooperativeDestroy(dCooperativeID cooperative);
dResourceRequirementsID dResourceRequirementsCreate(dCooperativeID cooperative);
void dResourceRequirementsDestroy(dResourceRequirementsID requirements);
dResourceRequirementsID dResourceRequirementsClone(dResourceRequirementsID requirements);
void dResourceRequirementsMergeIn(dResourceRequirementsID summaryRequirements, dResourceRequirementsID extraRequirements);
dResourceContainerID dResourceContainerAcquire(dResourceRequirementsID requirements);
void dResourceContainerDestroy(dResourceContainerID resources);
void dEstimateCooperativelyFactorLDLTResourceRequirements(dResourceRequirementsID requirements, unsigned maximalAllowedThreadCount, unsigned maximalRowCount);
void dCooperativelyFactorLDLT(dResourceContainerID resources, unsigned allowedThreadCount, dReal *A, dReal *d, unsigned rowCount, unsigned rowSkip);
void dEstimateCooperativelySolveLDLTResourceRequirements(dResourceRequirementsID requirements, unsigned maximalAllowedThreadCount, unsigned maximalRowCount);
void dCooperativelySolveLDLT(dResourceContainerID resources, unsigned allowedThreadCount, const dReal *L, const dReal *d, dReal *b, unsigned rowCount, unsigned rowSkip);
void dEstimateCooperativelySolveL1StraightResourceRequirements(dResourceRequirementsID requirements, unsigned maximalAllowedThreadCount, unsigned maximalRowCount);
void dCooperativelySolveL1Straight(dResourceContainerID resources, unsigned allowedThreadCount, const dReal *L, dReal *b, unsigned rowCount, unsigned rowSkip);
void dEstimateCooperativelySolveL1TransposedResourceRequirements(dResourceRequirementsID requirements, unsigned maximalAllowedThreadCount, unsigned maximalRowCount);
void dCooperativelySolveL1Transposed(dResourceContainerID resources, unsigned allowedThreadCount, const dReal *L, dReal *b, unsigned rowCount, unsigned rowSkip);
void dEstimateCooperativelyScaleVectorResourceRequirements(dResourceRequirementsID requirements, unsigned maximalAllowedThreadCount, unsigned maximalElementCount);
void dCooperativelyScaleVector(dResourceContainerID resources, unsigned allowedThreadCount, dReal *dataVector, const dReal *scaleVector, unsigned elementCount);
typedef struct dStopwatch {
  double time; 
  unsigned long cc[2]; 
} dStopwatch;
void dStopwatchReset (dStopwatch *);
void dStopwatchStart (dStopwatch *);
void dStopwatchStop (dStopwatch *);
double dStopwatchTime (dStopwatch *); 
void dTimerStart (const char *description); 
void dTimerNow (const char *description); 
void dTimerEnd(void);
typedef struct FILE FILE;
void dTimerReport (FILE *fout, int average);
double dTimerTicksPerSecond(void);
double dTimerResolution(void);
void dRSetIdentity (dMatrix3 R);
void dRFromAxisAndAngle (dMatrix3 R, dReal ax, dReal ay, dReal az, dReal angle);
void dRFromEulerAngles (dMatrix3 R, dReal phi, dReal theta, dReal psi);
void dRFrom2Axes (dMatrix3 R, dReal ax, dReal ay, dReal az, dReal bx, dReal by, dReal bz);
void dRFromZAxis (dMatrix3 R, dReal ax, dReal ay, dReal az);
void dQSetIdentity (dQuaternion q);
void dQFromAxisAndAngle (dQuaternion q, dReal ax, dReal ay, dReal az, dReal angle);
void dQMultiply0 (dQuaternion qa, const dQuaternion qb, const dQuaternion qc);
void dQMultiply1 (dQuaternion qa, const dQuaternion qb, const dQuaternion qc);
void dQMultiply2 (dQuaternion qa, const dQuaternion qb, const dQuaternion qc);
void dQMultiply3 (dQuaternion qa, const dQuaternion qb, const dQuaternion qc);
void dRfromQ (dMatrix3 R, const dQuaternion q);
void dQfromR (dQuaternion q, const dMatrix3 R);
void dDQfromW (dReal dq[4], const dVector3 w, const dQuaternion q);
struct dMass;
typedef struct dMass dMass;
int dMassCheck(const dMass *m);
void dMassSetZero (dMass *);
void dMassSetParameters (dMass *, dReal themass, dReal cgx, dReal cgy, dReal cgz, dReal I11, dReal I22, dReal I33, dReal I12, dReal I13, dReal I23);
void dMassSetSphere (dMass *, dReal density, dReal radius);
void dMassSetSphereTotal (dMass *, dReal total_mass, dReal radius);
void dMassSetCapsule (dMass *, dReal density, int direction, dReal radius, dReal length);
void dMassSetCapsuleTotal (dMass *, dReal total_mass, int direction, dReal radius, dReal length);
void dMassSetCylinder (dMass *, dReal density, int direction, dReal radius, dReal length);
void dMassSetCylinderTotal (dMass *, dReal total_mass, int direction, dReal radius, dReal length);
void dMassSetBox (dMass *, dReal density, dReal lx, dReal ly, dReal lz);
void dMassSetBoxTotal (dMass *, dReal total_mass, dReal lx, dReal ly, dReal lz);
void dMassSetTrimesh (dMass *, dReal density, dGeomID g);
void dMassSetTrimeshTotal (dMass *m, dReal total_mass, dGeomID g);
void dMassAdjust (dMass *, dReal newmass);
void dMassTranslate (dMass *, dReal x, dReal y, dReal z);
void dMassRotate (dMass *, const dMatrix3 R);
void dMassAdd (dMass *a, const dMass *b);
void dMassSetCappedCylinder(dMass *a, dReal b, int c, dReal d, dReal e);
void dMassSetCappedCylinderTotal(dMass *a, dReal b, int c, dReal d, dReal e);
struct dMass {
  dReal mass;
  dVector3 c;
  dMatrix3 I;
};
int dTestRand(void);
unsigned long dRand(void);
unsigned long dRandGetSeed(void);
void dRandSetSeed (unsigned long s);
int dRandInt (int n);
dReal dRandReal(void);
void dPrintMatrix (const dReal *A, int n, int m, const char *fmt, FILE *f);
void dMakeRandomVector (dReal *A, int n, dReal range);
void dMakeRandomMatrix (dReal *A, int n, int m, dReal range);
void dClearUpperTriangle (dReal *A, int n);
dReal dMaxDifference (const dReal *A, const dReal *B, int n, int m);
dReal dMaxDifferenceLowerTriangle (const dReal *A, const dReal *B, int n);
dWorldID dWorldCreate(void);
void dWorldDestroy (dWorldID world);
void dWorldSetData (dWorldID world, void* data);
void* dWorldGetData (dWorldID world);
void dWorldSetGravity (dWorldID, dReal x, dReal y, dReal z);
void dWorldGetGravity (dWorldID, dVector3 gravity);
void dWorldSetERP (dWorldID, dReal erp);
dReal dWorldGetERP (dWorldID);
void dWorldSetCFM (dWorldID, dReal cfm);
dReal dWorldGetCFM (dWorldID);
void dWorldSetStepIslandsProcessingMaxThreadCount(dWorldID w, unsigned count);
unsigned dWorldGetStepIslandsProcessingMaxThreadCount(dWorldID w);
int dWorldUseSharedWorkingMemory(dWorldID w, dWorldID from_world);
void dWorldCleanupWorkingMemory(dWorldID w);
typedef struct
{
  unsigned struct_size;
  float reserve_factor; 
  unsigned reserve_minimum;
} dWorldStepReserveInfo;
int dWorldSetStepMemoryReservationPolicy(dWorldID w, const dWorldStepReserveInfo *policyinfo);
typedef struct
{
  unsigned struct_size;
  void *(*alloc_block)(dsizeint block_size);
  void *(*shrink_block)(void *block_pointer, dsizeint block_current_size, dsizeint block_smaller_size);
  void (*free_block)(void *block_pointer, dsizeint block_current_size);
} dWorldStepMemoryFunctionsInfo;
int dWorldSetStepMemoryManager(dWorldID w, const dWorldStepMemoryFunctionsInfo *memfuncs);
void dWorldSetStepThreadingImplementation(dWorldID w, const dThreadingFunctionsInfo *functions_info, dThreadingImplementationID threading_impl);
int dWorldStep (dWorldID w, dReal stepsize);
int dWorldQuickStep (dWorldID w, dReal stepsize);
void dWorldImpulseToForce
(
dWorldID, dReal stepsize,
dReal ix, dReal iy, dReal iz, dVector3 force
 );
void dWorldSetQuickStepNumIterations (dWorldID, int num);
int dWorldGetQuickStepNumIterations (dWorldID);
void dWorldSetQuickStepW (dWorldID, dReal over_relaxation);
dReal dWorldGetQuickStepW (dWorldID);
void dWorldSetContactMaxCorrectingVel (dWorldID, dReal vel);
dReal dWorldGetContactMaxCorrectingVel (dWorldID);
void dWorldSetContactSurfaceLayer (dWorldID, dReal depth);
dReal dWorldGetContactSurfaceLayer (dWorldID);
dReal dWorldGetAutoDisableLinearThreshold (dWorldID);
void dWorldSetAutoDisableLinearThreshold (dWorldID, dReal linear_average_threshold);
dReal dWorldGetAutoDisableAngularThreshold (dWorldID);
void dWorldSetAutoDisableAngularThreshold (dWorldID, dReal angular_average_threshold);
int dWorldGetAutoDisableAverageSamplesCount (dWorldID);
void dWorldSetAutoDisableAverageSamplesCount (dWorldID, unsigned int average_samples_count );
int dWorldGetAutoDisableSteps (dWorldID);
void dWorldSetAutoDisableSteps (dWorldID, int steps);
dReal dWorldGetAutoDisableTime (dWorldID);
void dWorldSetAutoDisableTime (dWorldID, dReal time);
int dWorldGetAutoDisableFlag (dWorldID);
void dWorldSetAutoDisableFlag (dWorldID, int do_auto_disable);
dReal dWorldGetLinearDampingThreshold (dWorldID w);
void dWorldSetLinearDampingThreshold(dWorldID w, dReal threshold);
dReal dWorldGetAngularDampingThreshold (dWorldID w);
void dWorldSetAngularDampingThreshold(dWorldID w, dReal threshold);
dReal dWorldGetLinearDamping (dWorldID w);
void dWorldSetLinearDamping (dWorldID w, dReal scale);
dReal dWorldGetAngularDamping (dWorldID w);
void dWorldSetAngularDamping(dWorldID w, dReal scale);
void dWorldSetDamping(dWorldID w, dReal linear_scale, dReal angular_scale);
dReal dWorldGetMaxAngularSpeed (dWorldID w);
void dWorldSetMaxAngularSpeed (dWorldID w, dReal max_speed);
dReal dBodyGetAutoDisableLinearThreshold (dBodyID);
void dBodySetAutoDisableLinearThreshold (dBodyID, dReal linear_average_threshold);
dReal dBodyGetAutoDisableAngularThreshold (dBodyID);
void dBodySetAutoDisableAngularThreshold (dBodyID, dReal angular_average_threshold);
int dBodyGetAutoDisableAverageSamplesCount (dBodyID);
void dBodySetAutoDisableAverageSamplesCount (dBodyID, unsigned int average_samples_count);
int dBodyGetAutoDisableSteps (dBodyID);
void dBodySetAutoDisableSteps (dBodyID, int steps);
dReal dBodyGetAutoDisableTime (dBodyID);
void dBodySetAutoDisableTime (dBodyID, dReal time);
int dBodyGetAutoDisableFlag (dBodyID);
void dBodySetAutoDisableFlag (dBodyID, int do_auto_disable);
void dBodySetAutoDisableDefaults (dBodyID);
dWorldID dBodyGetWorld (dBodyID);
dBodyID dBodyCreate (dWorldID);
void dBodyDestroy (dBodyID);
void dBodySetData (dBodyID, void *data);
void *dBodyGetData (dBodyID);
void dBodySetPosition (dBodyID, dReal x, dReal y, dReal z);
void dBodySetRotation (dBodyID, const dMatrix3 R);
void dBodySetQuaternion (dBodyID, const dQuaternion q);
void dBodySetLinearVel (dBodyID, dReal x, dReal y, dReal z);
void dBodySetAngularVel (dBodyID, dReal x, dReal y, dReal z);
const dReal * dBodyGetPosition (dBodyID);
void dBodyCopyPosition (dBodyID body, dVector3 pos);
const dReal * dBodyGetRotation (dBodyID);
void dBodyCopyRotation (dBodyID, dMatrix3 R);
const dReal * dBodyGetQuaternion (dBodyID);
void dBodyCopyQuaternion(dBodyID body, dQuaternion quat);
const dReal * dBodyGetLinearVel (dBodyID);
const dReal * dBodyGetAngularVel (dBodyID);
void dBodySetMass (dBodyID, const dMass *mass);
void dBodyGetMass (dBodyID, dMass *mass);
void dBodyAddForce (dBodyID, dReal fx, dReal fy, dReal fz);
void dBodyAddTorque (dBodyID, dReal fx, dReal fy, dReal fz);
void dBodyAddRelForce (dBodyID, dReal fx, dReal fy, dReal fz);
void dBodyAddRelTorque (dBodyID, dReal fx, dReal fy, dReal fz);
void dBodyAddForceAtPos (dBodyID, dReal fx, dReal fy, dReal fz, dReal px, dReal py, dReal pz);
void dBodyAddForceAtRelPos (dBodyID, dReal fx, dReal fy, dReal fz, dReal px, dReal py, dReal pz);
void dBodyAddRelForceAtPos (dBodyID, dReal fx, dReal fy, dReal fz, dReal px, dReal py, dReal pz);
void dBodyAddRelForceAtRelPos (dBodyID, dReal fx, dReal fy, dReal fz, dReal px, dReal py, dReal pz);
const dReal * dBodyGetForce (dBodyID);
const dReal * dBodyGetTorque (dBodyID);
void dBodySetForce (dBodyID b, dReal x, dReal y, dReal z);
void dBodySetTorque (dBodyID b, dReal x, dReal y, dReal z);
void dBodyGetRelPointPos(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodyGetRelPointVel(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodyGetPointVel(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodyGetPosRelPoint(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodyVectorToWorld(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodyVectorFromWorld(dBodyID, dReal px, dReal py, dReal pz, dVector3 result);
void dBodySetFiniteRotationMode (dBodyID, int mode);
void dBodySetFiniteRotationAxis (dBodyID, dReal x, dReal y, dReal z);
int dBodyGetFiniteRotationMode (dBodyID);
void dBodyGetFiniteRotationAxis (dBodyID, dVector3 result);
int dBodyGetNumJoints (dBodyID b);
dJointID dBodyGetJoint (dBodyID, int index);
void dBodySetDynamic (dBodyID);
void dBodySetKinematic (dBodyID);
int dBodyIsKinematic (dBodyID);
void dBodyEnable (dBodyID);
void dBodyDisable (dBodyID);
int dBodyIsEnabled (dBodyID);
void dBodySetGravityMode (dBodyID b, int mode);
int dBodyGetGravityMode (dBodyID b);
void dBodySetMovedCallback(dBodyID b, void (*callback)(dBodyID));
dGeomID dBodyGetFirstGeom (dBodyID b);
dGeomID dBodyGetNextGeom (dGeomID g);
void dBodySetDampingDefaults(dBodyID b);
dReal dBodyGetLinearDamping (dBodyID b);
void dBodySetLinearDamping(dBodyID b, dReal scale);
dReal dBodyGetAngularDamping (dBodyID b);
void dBodySetAngularDamping(dBodyID b, dReal scale);
void dBodySetDamping(dBodyID b, dReal linear_scale, dReal angular_scale);
dReal dBodyGetLinearDampingThreshold (dBodyID b);
void dBodySetLinearDampingThreshold(dBodyID b, dReal threshold);
dReal dBodyGetAngularDampingThreshold (dBodyID b);
void dBodySetAngularDampingThreshold(dBodyID b, dReal threshold);
dReal dBodyGetMaxAngularSpeed (dBodyID b);
void dBodySetMaxAngularSpeed(dBodyID b, dReal max_speed);
int dBodyGetGyroscopicMode(dBodyID b);
void dBodySetGyroscopicMode(dBodyID b, int enabled);
dJointID dJointCreateBall (dWorldID, dJointGroupID);
dJointID dJointCreateHinge (dWorldID, dJointGroupID);
dJointID dJointCreateSlider (dWorldID, dJointGroupID);
dJointID dJointCreateContact (dWorldID, dJointGroupID, const dContact *);
dJointID dJointCreateHinge2 (dWorldID, dJointGroupID);
dJointID dJointCreateUniversal (dWorldID, dJointGroupID);
dJointID dJointCreatePR (dWorldID, dJointGroupID);
dJointID dJointCreatePU (dWorldID, dJointGroupID);
dJointID dJointCreatePiston (dWorldID, dJointGroupID);
dJointID dJointCreateFixed (dWorldID, dJointGroupID);
dJointID dJointCreateNull (dWorldID, dJointGroupID);
dJointID dJointCreateAMotor (dWorldID, dJointGroupID);
dJointID dJointCreateLMotor (dWorldID, dJointGroupID);
dJointID dJointCreatePlane2D (dWorldID, dJointGroupID);
dJointID dJointCreateDBall (dWorldID, dJointGroupID);
dJointID dJointCreateDHinge (dWorldID, dJointGroupID);
dJointID dJointCreateTransmission (dWorldID, dJointGroupID);
void dJointDestroy (dJointID);
dJointGroupID dJointGroupCreate (int max_size);
void dJointGroupDestroy (dJointGroupID);
void dJointGroupEmpty (dJointGroupID);
int dJointGetNumBodies(dJointID);
void dJointAttach (dJointID, dBodyID body1, dBodyID body2);
void dJointEnable (dJointID);
void dJointDisable (dJointID);
int dJointIsEnabled (dJointID);
void dJointSetData (dJointID, void *data);
void *dJointGetData (dJointID);
dJointType dJointGetType (dJointID);
dBodyID dJointGetBody (dJointID, int index);
void dJointSetFeedback (dJointID, dJointFeedback *);
dJointFeedback *dJointGetFeedback (dJointID);
void dJointSetBallAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetBallAnchor2 (dJointID, dReal x, dReal y, dReal z);
void dJointSetBallParam (dJointID, int parameter, dReal value);
void dJointSetHingeAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetHingeAnchorDelta (dJointID, dReal x, dReal y, dReal z, dReal ax, dReal ay, dReal az);
void dJointSetHingeAxis (dJointID, dReal x, dReal y, dReal z);
void dJointSetHingeAxisOffset (dJointID j, dReal x, dReal y, dReal z, dReal angle);
void dJointSetHingeParam (dJointID, int parameter, dReal value);
void dJointAddHingeTorque(dJointID joint, dReal torque);
void dJointSetSliderAxis (dJointID, dReal x, dReal y, dReal z);
void dJointSetSliderAxisDelta (dJointID, dReal x, dReal y, dReal z, dReal ax, dReal ay, dReal az);
void dJointSetSliderParam (dJointID, int parameter, dReal value);
void dJointAddSliderForce(dJointID joint, dReal force);
void dJointSetHinge2Anchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetHinge2Axes (dJointID j, const dReal *axis1, const dReal *axis2);
void dJointSetHinge2Axis1 (dJointID j, dReal x, dReal y, dReal z);
void dJointSetHinge2Axis2 (dJointID j, dReal x, dReal y, dReal z);
void dJointSetHinge2Param (dJointID, int parameter, dReal value);
void dJointAddHinge2Torques(dJointID joint, dReal torque1, dReal torque2);
void dJointSetUniversalAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetUniversalAxis1 (dJointID, dReal x, dReal y, dReal z);
void dJointSetUniversalAxis1Offset (dJointID, dReal x, dReal y, dReal z, dReal offset1, dReal offset2);
void dJointSetUniversalAxis2 (dJointID, dReal x, dReal y, dReal z);
void dJointSetUniversalAxis2Offset (dJointID, dReal x, dReal y, dReal z, dReal offset1, dReal offset2);
void dJointSetUniversalParam (dJointID, int parameter, dReal value);
void dJointAddUniversalTorques(dJointID joint, dReal torque1, dReal torque2);
void dJointSetPRAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetPRAxis1 (dJointID, dReal x, dReal y, dReal z);
void dJointSetPRAxis2 (dJointID, dReal x, dReal y, dReal z);
void dJointSetPRParam (dJointID, int parameter, dReal value);
void dJointAddPRTorque (dJointID j, dReal torque);
void dJointSetPUAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetPUAnchorDelta (dJointID, dReal x, dReal y, dReal z, dReal dx, dReal dy, dReal dz);
void dJointSetPUAnchorOffset (dJointID, dReal x, dReal y, dReal z, dReal dx, dReal dy, dReal dz);
void dJointSetPUAxis1 (dJointID, dReal x, dReal y, dReal z);
void dJointSetPUAxis2 (dJointID, dReal x, dReal y, dReal z);
void dJointSetPUAxis3 (dJointID, dReal x, dReal y, dReal z);
void dJointSetPUAxisP (dJointID id, dReal x, dReal y, dReal z);
void dJointSetPUParam (dJointID, int parameter, dReal value);
void dJointSetPistonAnchor (dJointID, dReal x, dReal y, dReal z);
void dJointSetPistonAnchorOffset(dJointID j, dReal x, dReal y, dReal z, dReal dx, dReal dy, dReal dz);
void dJointSetPistonAxis (dJointID, dReal x, dReal y, dReal z);
void dJointSetPistonAxisDelta (dJointID j, dReal x, dReal y, dReal z, dReal ax, dReal ay, dReal az);
void dJointSetPistonParam (dJointID, int parameter, dReal value);
void dJointAddPistonForce (dJointID joint, dReal force);
void dJointSetFixed (dJointID);
void dJointSetFixedParam (dJointID, int parameter, dReal value);
void dJointSetAMotorNumAxes (dJointID, int num);
void dJointSetAMotorAxis (dJointID, int anum, int rel, dReal x, dReal y, dReal z);
void dJointSetAMotorAngle (dJointID, int anum, dReal angle);
void dJointSetAMotorParam (dJointID, int parameter, dReal value);
void dJointSetAMotorMode (dJointID, int mode);
void dJointAddAMotorTorques (dJointID, dReal torque1, dReal torque2, dReal torque3);
void dJointSetLMotorNumAxes (dJointID, int num);
void dJointSetLMotorAxis (dJointID, int anum, int rel, dReal x, dReal y, dReal z);
void dJointSetLMotorParam (dJointID, int parameter, dReal value);
void dJointSetPlane2DXParam (dJointID, int parameter, dReal value);
void dJointSetPlane2DYParam (dJointID, int parameter, dReal value);
void dJointSetPlane2DAngleParam (dJointID, int parameter, dReal value);
void dJointGetBallAnchor (dJointID, dVector3 result);
void dJointGetBallAnchor2 (dJointID, dVector3 result);
dReal dJointGetBallParam (dJointID, int parameter);
void dJointGetHingeAnchor (dJointID, dVector3 result);
void dJointGetHingeAnchor2 (dJointID, dVector3 result);
void dJointGetHingeAxis (dJointID, dVector3 result);
dReal dJointGetHingeParam (dJointID, int parameter);
dReal dJointGetHingeAngle (dJointID);
dReal dJointGetHingeAngleRate (dJointID);
dReal dJointGetSliderPosition (dJointID);
dReal dJointGetSliderPositionRate (dJointID);
void dJointGetSliderAxis (dJointID, dVector3 result);
dReal dJointGetSliderParam (dJointID, int parameter);
void dJointGetHinge2Anchor (dJointID, dVector3 result);
void dJointGetHinge2Anchor2 (dJointID, dVector3 result);
void dJointGetHinge2Axis1 (dJointID, dVector3 result);
void dJointGetHinge2Axis2 (dJointID, dVector3 result);
dReal dJointGetHinge2Param (dJointID, int parameter);
dReal dJointGetHinge2Angle1 (dJointID);
dReal dJointGetHinge2Angle2 (dJointID);
dReal dJointGetHinge2Angle1Rate (dJointID);
dReal dJointGetHinge2Angle2Rate (dJointID);
void dJointGetUniversalAnchor (dJointID, dVector3 result);
void dJointGetUniversalAnchor2 (dJointID, dVector3 result);
void dJointGetUniversalAxis1 (dJointID, dVector3 result);
void dJointGetUniversalAxis2 (dJointID, dVector3 result);
dReal dJointGetUniversalParam (dJointID, int parameter);
void dJointGetUniversalAngles (dJointID, dReal *angle1, dReal *angle2);
dReal dJointGetUniversalAngle1 (dJointID);
dReal dJointGetUniversalAngle2 (dJointID);
dReal dJointGetUniversalAngle1Rate (dJointID);
dReal dJointGetUniversalAngle2Rate (dJointID);
void dJointGetPRAnchor (dJointID, dVector3 result);
dReal dJointGetPRPosition (dJointID);
dReal dJointGetPRPositionRate (dJointID);
dReal dJointGetPRAngle (dJointID);
dReal dJointGetPRAngleRate (dJointID);
void dJointGetPRAxis1 (dJointID, dVector3 result);
void dJointGetPRAxis2 (dJointID, dVector3 result);
dReal dJointGetPRParam (dJointID, int parameter);
void dJointGetPUAnchor (dJointID, dVector3 result);
dReal dJointGetPUPosition (dJointID);
dReal dJointGetPUPositionRate (dJointID);
void dJointGetPUAxis1 (dJointID, dVector3 result);
void dJointGetPUAxis2 (dJointID, dVector3 result);
void dJointGetPUAxis3 (dJointID, dVector3 result);
void dJointGetPUAxisP (dJointID id, dVector3 result);
void dJointGetPUAngles (dJointID, dReal *angle1, dReal *angle2);
dReal dJointGetPUAngle1 (dJointID);
dReal dJointGetPUAngle1Rate (dJointID);
dReal dJointGetPUAngle2 (dJointID);
dReal dJointGetPUAngle2Rate (dJointID);
dReal dJointGetPUParam (dJointID, int parameter);
dReal dJointGetPistonPosition (dJointID);
dReal dJointGetPistonPositionRate (dJointID);
dReal dJointGetPistonAngle (dJointID);
dReal dJointGetPistonAngleRate (dJointID);
void dJointGetPistonAnchor (dJointID, dVector3 result);
void dJointGetPistonAnchor2 (dJointID, dVector3 result);
void dJointGetPistonAxis (dJointID, dVector3 result);
dReal dJointGetPistonParam (dJointID, int parameter);
int dJointGetAMotorNumAxes (dJointID);
void dJointGetAMotorAxis (dJointID, int anum, dVector3 result);
int dJointGetAMotorAxisRel (dJointID, int anum);
dReal dJointGetAMotorAngle (dJointID, int anum);
dReal dJointGetAMotorAngleRate (dJointID, int anum);
dReal dJointGetAMotorParam (dJointID, int parameter);
int dJointGetAMotorMode (dJointID);
int dJointGetLMotorNumAxes (dJointID);
void dJointGetLMotorAxis (dJointID, int anum, dVector3 result);
dReal dJointGetLMotorParam (dJointID, int parameter);
dReal dJointGetFixedParam (dJointID, int parameter);
void dJointGetTransmissionContactPoint1(dJointID, dVector3 result);
void dJointGetTransmissionContactPoint2(dJointID, dVector3 result);
void dJointSetTransmissionAxis1(dJointID, dReal x, dReal y, dReal z);
void dJointGetTransmissionAxis1(dJointID, dVector3 result);
void dJointSetTransmissionAxis2(dJointID, dReal x, dReal y, dReal z);
void dJointGetTransmissionAxis2(dJointID, dVector3 result);
void dJointSetTransmissionAnchor1(dJointID, dReal x, dReal y, dReal z);
void dJointGetTransmissionAnchor1(dJointID, dVector3 result);
void dJointSetTransmissionAnchor2(dJointID, dReal x, dReal y, dReal z);
void dJointGetTransmissionAnchor2(dJointID, dVector3 result);
void dJointSetTransmissionParam(dJointID, int parameter, dReal value);
dReal dJointGetTransmissionParam(dJointID, int parameter);
void dJointSetTransmissionMode( dJointID j, int mode );
int dJointGetTransmissionMode( dJointID j );
void dJointSetTransmissionRatio( dJointID j, dReal ratio );
dReal dJointGetTransmissionRatio( dJointID j );
void dJointSetTransmissionAxis( dJointID j, dReal x, dReal y, dReal z );
void dJointGetTransmissionAxis( dJointID j, dVector3 result );
dReal dJointGetTransmissionAngle1( dJointID j );
dReal dJointGetTransmissionAngle2( dJointID j );
dReal dJointGetTransmissionRadius1( dJointID j );
dReal dJointGetTransmissionRadius2( dJointID j );
void dJointSetTransmissionRadius1( dJointID j, dReal radius );
void dJointSetTransmissionRadius2( dJointID j, dReal radius );
dReal dJointGetTransmissionBacklash( dJointID j );
void dJointSetTransmissionBacklash( dJointID j, dReal backlash );
void dJointSetDBallAnchor1(dJointID, dReal x, dReal y, dReal z);
void dJointSetDBallAnchor2(dJointID, dReal x, dReal y, dReal z);
void dJointGetDBallAnchor1(dJointID, dVector3 result);
void dJointGetDBallAnchor2(dJointID, dVector3 result);
dReal dJointGetDBallDistance(dJointID);
void dJointSetDBallDistance(dJointID, dReal dist);
void dJointSetDBallParam(dJointID, int parameter, dReal value);
dReal dJointGetDBallParam(dJointID, int parameter);
void dJointSetDHingeAxis(dJointID, dReal x, dReal y, dReal z);
void dJointGetDHingeAxis(dJointID, dVector3 result);
void dJointSetDHingeAnchor1(dJointID, dReal x, dReal y, dReal z);
void dJointSetDHingeAnchor2(dJointID, dReal x, dReal y, dReal z);
void dJointGetDHingeAnchor1(dJointID, dVector3 result);
void dJointGetDHingeAnchor2(dJointID, dVector3 result);
dReal dJointGetDHingeDistance(dJointID);
void dJointSetDHingeParam(dJointID, int parameter, dReal value);
dReal dJointGetDHingeParam(dJointID, int parameter);
dJointID dConnectingJoint (dBodyID, dBodyID);
int dConnectingJointList (dBodyID, dBodyID, dJointID*);
int dAreConnected (dBodyID, dBodyID);
int dAreConnectedExcluding (dBodyID body1, dBodyID body2, int joint_type);
struct dContactGeom;
typedef void dNearCallback (void *data, dGeomID o1, dGeomID o2);
dSpaceID dSimpleSpaceCreate (dSpaceID space);
dSpaceID dHashSpaceCreate (dSpaceID space);
dSpaceID dQuadTreeSpaceCreate (dSpaceID space, const dVector3 Center, const dVector3 Extents, int Depth);
dSpaceID dSweepAndPruneSpaceCreate( dSpaceID space, int axisorder );
void dSpaceDestroy (dSpaceID);
void dHashSpaceSetLevels (dSpaceID space, int minlevel, int maxlevel);
void dHashSpaceGetLevels (dSpaceID space, int *minlevel, int *maxlevel);
void dSpaceSetCleanup (dSpaceID space, int mode);
int dSpaceGetCleanup (dSpaceID space);
void dSpaceSetSublevel (dSpaceID space, int sublevel);
int dSpaceGetSublevel (dSpaceID space);
void dSpaceSetManualCleanup (dSpaceID space, int mode);
int dSpaceGetManualCleanup (dSpaceID space);
void dSpaceAdd (dSpaceID, dGeomID);
void dSpaceRemove (dSpaceID, dGeomID);
int dSpaceQuery (dSpaceID, dGeomID);
void dSpaceClean (dSpaceID);
int dSpaceGetNumGeoms (dSpaceID);
dGeomID dSpaceGetGeom (dSpaceID, int i);
int dSpaceGetClass(dSpaceID space);
void dGeomDestroy (dGeomID geom);
void dGeomSetData (dGeomID geom, void* data);
void *dGeomGetData (dGeomID geom);
void dGeomSetBody (dGeomID geom, dBodyID body);
dBodyID dGeomGetBody (dGeomID geom);
void dGeomSetPosition (dGeomID geom, dReal x, dReal y, dReal z);
void dGeomSetRotation (dGeomID geom, const dMatrix3 R);
void dGeomSetQuaternion (dGeomID geom, const dQuaternion Q);
const dReal * dGeomGetPosition (dGeomID geom);
void dGeomCopyPosition (dGeomID geom, dVector3 pos);
const dReal * dGeomGetRotation (dGeomID geom);
void dGeomCopyRotation(dGeomID geom, dMatrix3 R);
void dGeomGetQuaternion (dGeomID geom, dQuaternion result);
void dGeomGetAABB (dGeomID geom, dReal aabb[6]);
int dGeomIsSpace (dGeomID geom);
dSpaceID dGeomGetSpace (dGeomID);
int dGeomGetClass (dGeomID geom);
void dGeomSetCategoryBits (dGeomID geom, unsigned long bits);
void dGeomSetCollideBits (dGeomID geom, unsigned long bits);
unsigned long dGeomGetCategoryBits (dGeomID);
unsigned long dGeomGetCollideBits (dGeomID);
void dGeomEnable (dGeomID geom);
void dGeomDisable (dGeomID geom);
int dGeomIsEnabled (dGeomID geom);
enum
{
dGeomCommonControlClass = 0,
dGeomColliderControlClass = 1
};
enum
{
dGeomCommonAnyControlCode = 0,
dGeomColliderSetMergeSphereContactsControlCode = 1,
dGeomColliderGetMergeSphereContactsControlCode = 2
};
enum
{
    dGeomColliderMergeContactsValue__Default = 0, 
dGeomColliderMergeContactsValue_None = 1,
dGeomColliderMergeContactsValue_Normals = 2,
dGeomColliderMergeContactsValue_Full = 3
};
int dGeomLowLevelControl (dGeomID geom, int controlClass, int controlCode, void *dataValue, int *dataSize);
void dGeomGetRelPointPos(dGeomID geom, dReal px, dReal py, dReal pz, dVector3 result);
void dGeomGetPosRelPoint(dGeomID geom, dReal px, dReal py, dReal pz, dVector3 result);
void dGeomVectorToWorld(dGeomID geom, dReal px, dReal py, dReal pz, dVector3 result);
void dGeomVectorFromWorld(dGeomID geom, dReal px, dReal py, dReal pz, dVector3 result);
void dGeomSetOffsetPosition (dGeomID geom, dReal x, dReal y, dReal z);
void dGeomSetOffsetRotation (dGeomID geom, const dMatrix3 R);
void dGeomSetOffsetQuaternion (dGeomID geom, const dQuaternion Q);
void dGeomSetOffsetWorldPosition (dGeomID geom, dReal x, dReal y, dReal z);
void dGeomSetOffsetWorldRotation (dGeomID geom, const dMatrix3 R);
void dGeomSetOffsetWorldQuaternion (dGeomID geom, const dQuaternion);
void dGeomClearOffset(dGeomID geom);
int dGeomIsOffset(dGeomID geom);
const dReal * dGeomGetOffsetPosition (dGeomID geom);
void dGeomCopyOffsetPosition (dGeomID geom, dVector3 pos);
const dReal * dGeomGetOffsetRotation (dGeomID geom);
void dGeomCopyOffsetRotation (dGeomID geom, dMatrix3 R);
void dGeomGetOffsetQuaternion (dGeomID geom, dQuaternion result);
int dCollide (dGeomID o1, dGeomID o2, int flags, dContactGeom *contact, int skip);
void dSpaceCollide (dSpaceID space, void *data, dNearCallback *callback);
void dSpaceCollide2 (dGeomID space1, dGeomID space2, void *data, dNearCallback *callback);
enum {
  dMaxUserClasses = 4
};
enum {
  dSphereClass = 0, dBoxClass, dCapsuleClass, dCylinderClass, dPlaneClass, dRayClass, dConvexClass, dGeomTransformClass, dTriMeshClass, dHeightfieldClass, dFirstSpaceClass, dSimpleSpaceClass = dFirstSpaceClass, dHashSpaceClass, dSweepAndPruneSpaceClass, 
  dQuadTreeSpaceClass, dLastSpaceClass = dQuadTreeSpaceClass, dFirstUserClass, dLastUserClass = dFirstUserClass + dMaxUserClasses - 1, dGeomNumClasses
};
dGeomID dCreateSphere (dSpaceID space, dReal radius);
void dGeomSphereSetRadius (dGeomID sphere, dReal radius);
dReal dGeomSphereGetRadius (dGeomID sphere);
dReal dGeomSpherePointDepth (dGeomID sphere, dReal x, dReal y, dReal z);
dGeomID dCreateConvex (dSpaceID space, const dReal *_planes, unsigned int _planecount, const dReal *_points, unsigned int _pointcount, const unsigned int *_polygons);
void dGeomSetConvex (dGeomID g, const dReal *_planes, unsigned int _count, const dReal *_points, unsigned int _pointcount, const unsigned int *_polygons);
dGeomID dCreateBox (dSpaceID space, dReal lx, dReal ly, dReal lz);
void dGeomBoxSetLengths (dGeomID box, dReal lx, dReal ly, dReal lz);
void dGeomBoxGetLengths (dGeomID box, dVector3 result);
dReal dGeomBoxPointDepth (dGeomID box, dReal x, dReal y, dReal z);
dGeomID dCreatePlane (dSpaceID space, dReal a, dReal b, dReal c, dReal d);
void dGeomPlaneSetParams (dGeomID plane, dReal a, dReal b, dReal c, dReal d);
void dGeomPlaneGetParams (dGeomID plane, dVector4 result);
dReal dGeomPlanePointDepth (dGeomID plane, dReal x, dReal y, dReal z);
dGeomID dCreateCapsule (dSpaceID space, dReal radius, dReal length);
void dGeomCapsuleSetParams (dGeomID ccylinder, dReal radius, dReal length);
void dGeomCapsuleGetParams (dGeomID ccylinder, dReal *radius, dReal *length);
dReal dGeomCapsulePointDepth (dGeomID ccylinder, dReal x, dReal y, dReal z);
dGeomID dCreateCylinder (dSpaceID space, dReal radius, dReal length);
void dGeomCylinderSetParams (dGeomID cylinder, dReal radius, dReal length);
void dGeomCylinderGetParams (dGeomID cylinder, dReal *radius, dReal *length);
dGeomID dCreateRay (dSpaceID space, dReal length);
void dGeomRaySetLength (dGeomID ray, dReal length);
dReal dGeomRayGetLength (dGeomID ray);
void dGeomRaySet (dGeomID ray, dReal px, dReal py, dReal pz, dReal dx, dReal dy, dReal dz);
void dGeomRayGet (dGeomID ray, dVector3 start, dVector3 dir);
void dGeomRaySetParams (dGeomID g, int FirstContact, int BackfaceCull);
void dGeomRayGetParams (dGeomID g, int *FirstContact, int *BackfaceCull);
void dGeomRaySetFirstContact (dGeomID g, int firstContact);
int dGeomRayGetFirstContact (dGeomID g);
void dGeomRaySetBackfaceCull (dGeomID g, int backfaceCull);
int dGeomRayGetBackfaceCull (dGeomID g);
void dGeomRaySetClosestHit (dGeomID g, int closestHit);
int dGeomRayGetClosestHit (dGeomID g);
struct dxTriMeshData;
typedef struct dxTriMeshData* dTriMeshDataID;
typedef enum
{
    dMTV__MIN, dMTV_FIRST = dMTV__MIN, dMTV_SECOND, dMTV_THIRD, dMTV__MAX,
} dMeshTriangleVertex;
dTriMeshDataID dGeomTriMeshDataCreate(void);
void dGeomTriMeshDataDestroy(dTriMeshDataID g);
enum
{
    dTRIMESHDATA__MIN, dTRIMESHDATA_FACE_NORMALS = dTRIMESHDATA__MIN, dTRIMESHDATA_USE_FLAGS, dTRIMESHDATA__MAX, 
    TRIMESH_FACE_NORMALS = dTRIMESHDATA_FACE_NORMALS,
};
enum
{
    dMESHDATAUSE_EDGE1 = 0x01, dMESHDATAUSE_EDGE2 = 0x02, dMESHDATAUSE_EDGE3 = 0x04, dMESHDATAUSE_VERTEX1 = 0x08, dMESHDATAUSE_VERTEX2 = 0x10, dMESHDATAUSE_VERTEX3 = 0x20,
};
void dGeomTriMeshDataSet(dTriMeshDataID g, int data_id, void *in_data);
void *dGeomTriMeshDataGet2(dTriMeshDataID g, int data_id, dsizeint *pout_size);
void dGeomTriMeshSetLastTransform( dGeomID g, const dMatrix4 last_trans );
const dReal* dGeomTriMeshGetLastTransform( dGeomID g );
void dGeomTriMeshDataBuildSingle(dTriMeshDataID g, const void* Vertices, int VertexStride, int VertexCount, const void* Indices, int IndexCount, int TriStride);
void dGeomTriMeshDataBuildSingle1(dTriMeshDataID g, const void* Vertices, int VertexStride, int VertexCount, const void* Indices, int IndexCount, int TriStride, const void* Normals);
void dGeomTriMeshDataBuildDouble(dTriMeshDataID g, const void* Vertices, int VertexStride, int VertexCount, const void* Indices, int IndexCount, int TriStride);
void dGeomTriMeshDataBuildDouble1(dTriMeshDataID g, const void* Vertices, int VertexStride, int VertexCount, const void* Indices, int IndexCount, int TriStride, const void* Normals);
void dGeomTriMeshDataBuildSimple(dTriMeshDataID g, const dReal* Vertices, int VertexCount, const dTriIndex* Indices, int IndexCount);
void dGeomTriMeshDataBuildSimple1(dTriMeshDataID g, const dReal* Vertices, int VertexCount, const dTriIndex* Indices, int IndexCount, const int* Normals);
enum
{
    dTRIDATAPREPROCESS_BUILD__MIN, dTRIDATAPREPROCESS_BUILD_CONCAVE_EDGES = dTRIDATAPREPROCESS_BUILD__MIN, 
    dTRIDATAPREPROCESS_BUILD_FACE_ANGLES, 
    dTRIDATAPREPROCESS_BUILD__MAX,
};
enum
{
    dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA__MIN, dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA_BYTE_POSITIVE = dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA__MIN, 
    dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA_BYTE_ALL, 
    dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA_WORD_ALL, 
    dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA__MAX, dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA__DEFAULT = dTRIDATAPREPROCESS_FACE_ANGLES_EXTRA_BYTE_POSITIVE, 
};
int dGeomTriMeshDataPreprocess2(dTriMeshDataID g, unsigned int buildRequestFlags, const dintptr *requestExtraData);
int dGeomTriMeshDataPreprocess(dTriMeshDataID g);
void dGeomTriMeshDataGetBuffer(dTriMeshDataID g, unsigned char** buf, int* bufLen);
void dGeomTriMeshDataSetBuffer(dTriMeshDataID g, unsigned char* buf);
typedef int dTriCallback(dGeomID TriMesh, dGeomID RefObject, int TriangleIndex);
void dGeomTriMeshSetCallback(dGeomID g, dTriCallback* Callback);
dTriCallback* dGeomTriMeshGetCallback(dGeomID g);
typedef void dTriArrayCallback(dGeomID TriMesh, dGeomID RefObject, const int* TriIndices, int TriCount);
void dGeomTriMeshSetArrayCallback(dGeomID g, dTriArrayCallback* ArrayCallback);
dTriArrayCallback* dGeomTriMeshGetArrayCallback(dGeomID g);
typedef int dTriRayCallback(dGeomID TriMesh, dGeomID Ray, int TriangleIndex, dReal u, dReal v);
void dGeomTriMeshSetRayCallback(dGeomID g, dTriRayCallback* Callback);
dTriRayCallback* dGeomTriMeshGetRayCallback(dGeomID g);
typedef int dTriTriMergeCallback(dGeomID TriMesh, int FirstTriangleIndex, int SecondTriangleIndex);
void dGeomTriMeshSetTriMergeCallback(dGeomID g, dTriTriMergeCallback* Callback);
dTriTriMergeCallback* dGeomTriMeshGetTriMergeCallback(dGeomID g);
dGeomID dCreateTriMesh(dSpaceID space, dTriMeshDataID Data, dTriCallback* Callback, dTriArrayCallback* ArrayCallback, dTriRayCallback* RayCallback);
void dGeomTriMeshSetData(dGeomID g, dTriMeshDataID Data);
dTriMeshDataID dGeomTriMeshGetData(dGeomID g);
void dGeomTriMeshEnableTC(dGeomID g, int geomClass, int enable);
int dGeomTriMeshIsTCEnabled(dGeomID g, int geomClass);
void dGeomTriMeshClearTCCache(dGeomID g);
dTriMeshDataID dGeomTriMeshGetTriMeshDataID(dGeomID g);
void dGeomTriMeshGetTriangle(dGeomID g, int Index, dVector3* v0, dVector3* v1, dVector3* v2);
void dGeomTriMeshGetPoint(dGeomID g, int Index, dReal u, dReal v, dVector3 Out);
int dGeomTriMeshGetTriangleCount (dGeomID g);
void dGeomTriMeshDataUpdate(dTriMeshDataID g);
dGeomID dCreateGeomTransform (dSpaceID space);
void dGeomTransformSetGeom (dGeomID g, dGeomID obj);
dGeomID dGeomTransformGetGeom (dGeomID g);
void dGeomTransformSetCleanup (dGeomID g, int mode);
int dGeomTransformGetCleanup (dGeomID g);
void dGeomTransformSetInfo (dGeomID g, int mode);
int dGeomTransformGetInfo (dGeomID g);
struct dxHeightfieldData;
typedef struct dxHeightfieldData* dHeightfieldDataID;
typedef dReal dHeightfieldGetHeight( void* p_user_data, int x, int z );
dGeomID dCreateHeightfield( dSpaceID space, dHeightfieldDataID data, int bPlaceable );
dHeightfieldDataID dGeomHeightfieldDataCreate(void);
void dGeomHeightfieldDataDestroy( dHeightfieldDataID d );
void dGeomHeightfieldDataBuildCallback( dHeightfieldDataID d, void* pUserData, dHeightfieldGetHeight* pCallback, dReal width, dReal depth, int widthSamples, int depthSamples, dReal scale, dReal offset, dReal thickness, int bWrap );
void dGeomHeightfieldDataBuildByte( dHeightfieldDataID d, const unsigned char* pHeightData, int bCopyHeightData, dReal width, dReal depth, int widthSamples, int depthSamples, dReal scale, dReal offset, dReal thickness, int bWrap );
void dGeomHeightfieldDataBuildShort( dHeightfieldDataID d, const short* pHeightData, int bCopyHeightData, dReal width, dReal depth, int widthSamples, int depthSamples, dReal scale, dReal offset, dReal thickness, int bWrap );
void dGeomHeightfieldDataBuildSingle( dHeightfieldDataID d, const float* pHeightData, int bCopyHeightData, dReal width, dReal depth, int widthSamples, int depthSamples, dReal scale, dReal offset, dReal thickness, int bWrap );
void dGeomHeightfieldDataBuildDouble( dHeightfieldDataID d, const double* pHeightData, int bCopyHeightData, dReal width, dReal depth, int widthSamples, int depthSamples, dReal scale, dReal offset, dReal thickness, int bWrap );
void dGeomHeightfieldDataSetBounds( dHeightfieldDataID d, dReal minHeight, dReal maxHeight );
void dGeomHeightfieldSetHeightfieldData( dGeomID g, dHeightfieldDataID d );
dHeightfieldDataID dGeomHeightfieldGetHeightfieldData( dGeomID g );
void dClosestLineSegmentPoints (const dVector3 a1, const dVector3 a2, const dVector3 b1, const dVector3 b2, dVector3 cp1, dVector3 cp2);
int dBoxTouchesBox (const dVector3 _p1, const dMatrix3 R1, const dVector3 side1, const dVector3 _p2, const dMatrix3 R2, const dVector3 side2);
int dBoxBox (const dVector3 p1, const dMatrix3 R1, const dVector3 side1, const dVector3 p2, const dMatrix3 R2, const dVector3 side2, dVector3 normal, dReal *depth, int *return_code, int flags, dContactGeom *contact, int skip);
void dInfiniteAABB (dGeomID geom, dReal aabb[6]);
typedef void dGetAABBFn (dGeomID, dReal aabb[6]);
typedef int dColliderFn (dGeomID o1, dGeomID o2, int flags, dContactGeom *contact, int skip);
typedef dColliderFn * dGetColliderFnFn (int num);
typedef void dGeomDtorFn (dGeomID o);
typedef int dAABBTestFn (dGeomID o1, dGeomID o2, dReal aabb[6]);
typedef struct dGeomClass {
  int bytes;
  dGetColliderFnFn *collider;
  dGetAABBFn *aabb;
  dAABBTestFn *aabb_test;
  dGeomDtorFn *dtor;
} dGeomClass;
int dCreateGeomClass (const dGeomClass *classptr);
void * dGeomGetClassData (dGeomID);
dGeomID dCreateGeom (int classnum);
void dSetColliderOverride (int i, int j, dColliderFn *fn);
struct dxThreadingThreadPool;
typedef struct dxThreadingThreadPool *dThreadingThreadPoolID;
dThreadingImplementationID dThreadingAllocateSelfThreadedImplementation();
dThreadingImplementationID dThreadingAllocateMultiThreadedImplementation();
const dThreadingFunctionsInfo *dThreadingImplementationGetFunctions(dThreadingImplementationID impl);
void dThreadingImplementationShutdownProcessing(dThreadingImplementationID impl);
void dThreadingImplementationCleanupForRestart(dThreadingImplementationID impl);
void dThreadingFreeImplementation(dThreadingImplementationID impl);
typedef void (dThreadReadyToServeCallback)(void *callback_context);
void dExternalThreadingServeMultiThreadedImplementation(dThreadingImplementationID impl, dThreadReadyToServeCallback *readiness_callback, void *callback_context);
dThreadingThreadPoolID dThreadingAllocateThreadPool(unsigned thread_count, dsizeint stack_size, unsigned int ode_data_allocate_flags, void *reserved);
void dThreadingThreadPoolServeMultiThreadedImplementation(dThreadingThreadPoolID pool, dThreadingImplementationID impl);
void dThreadingThreadPoolWaitIdleState(dThreadingThreadPoolID pool);
void dThreadingFreeThreadPool(dThreadingThreadPoolID pool);
void dWorldExportDIF (dWorldID w, FILE *file, const char *world_name);
// Non-ODE functions
FILE *fopen(const char *filename, const char *mode);
int fclose(FILE *stream);
