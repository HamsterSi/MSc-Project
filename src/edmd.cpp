
#include "edmd.hpp"

/*******************************************************************************
 *  -/_|:|_|_\-
 *
 *  File:           qcprot.c
 *  Version:        1.4
 *
 *  Function:       Rapid calculation of the least-squares rotation using a
 *                  quaternion-based characteristic polynomial and
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 *
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for macromolecular
 *      superpositions."
 *      Journal of Computational Chemistry 31(7):1561-1563.
 *
 *
 *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *    of conditions and the following disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *    endorse or promote products derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
 *                    If trying all rows of the adjoint still gives too small
 *                    qsqr, then just return identity matrix. (DLT)
 *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
 *                    invalid mem access
 *    2011/02/21      Made CenterCoords use weights
 *    2011/05/02      Finally changed CenterCoords declaration in qcprot.h
 *                    Also changed some functions to static
 *    2011/07/08      put in fabs() to fix taking sqrt of small neg numbers, fp error
 *    2012/07/26      minor changes to comments and main.c, more info (v.1.4)
 *
 ******************************************************************************/

static double
InnerProduct(double *A, double **coords1, double **coords2, const int len, const double *weight)
{
    double          x1, x2, y1, y2, z1, z2;
    int             i;
    const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
    const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
    double          G1 = 0.0, G2 = 0.0;
    
    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;
    
    if (weight != NULL)
    {
        for (i = 0; i < len; ++i)
        {
            x1 = weight[i] * fx1[i];
            y1 = weight[i] * fy1[i];
            z1 = weight[i] * fz1[i];
            
            G1 += x1 * fx1[i] + y1 * fy1[i] + z1 * fz1[i];
            
            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];
            
            G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);
            
            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);
            
            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);
            
            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);
        }
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            x1 = fx1[i];
            y1 = fy1[i];
            z1 = fz1[i];
            
            G1 += x1 * x1 + y1 * y1 + z1 * z1;
            
            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];
            
            G2 += (x2 * x2 + y2 * y2 + z2 * z2);
            
            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);
            
            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);
            
            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);
        }
    }
    
    return (G1 + G2) * 0.5;
}


int
FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double C[4];
    int i;
    double mxEigenV;
    double oldg = 0.0;
    double b, a, delta, rms, qsqr;
    double q1, q2, q3, q4, normq;
    double a11, a12, a13, a14, a21, a22, a23, a24;
    double a31, a32, a33, a34, a41, a42, a43, a44;
    double a2, x2, y2, z2;
    double xy, az, zx, ay, yz, ax;
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132;
    double evecprec = 1e-6;
    double evalprec = 1e-11;
    
    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];
    
    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;
    
    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;
    
    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;
    
    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;
    
    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);
    
    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;
    
    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
    + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
    + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
    + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
    + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
    
    /* Newton-Raphson */
    mxEigenV = E0;
    for (i = 0; i < 50; ++i)
    {
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        /* printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV); */
        if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
            break;
    }
    
    if (i == 50)
        fprintf(stderr,"\nMore than %d iterations needed!\n", i);
    
    /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
    rms = sqrt(fabs(2.0 * (E0 - mxEigenV)/len));
    (*rmsd) = rms;
    /* printf("\n\n %16g %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); */
    
    if (minScore > 0)
        if (rms < minScore)
            return (-1); // Don't bother with rotation.
    
    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;
    
    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;
    
    /* The following code tries to calculate another column in the adjoint matrix when the norm of the
     current column is too small.
     Usually this block will never be activated.  To be absolutely safe this should be
     uncommented, but it is most likely unnecessary.
     */
    if (qsqr < evecprec)
    {
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;
        
        if (qsqr < evecprec)
        {
            double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
            double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
            double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;
            
            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;
            
            if (qsqr < evecprec)
            {
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr < evecprec)
                {
                    /* if qsqr is still too small, return the identity matrix. */
                    rot[0] = rot[4] = rot[8] = 1.0;
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;
                    
                    return(0);
                }
            }
        }
    }
    
    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;
    
    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;
    
    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;
    
    rot[0] = a2 + x2 - y2 - z2;
    rot[1] = 2 * (xy + az);
    rot[2] = 2 * (zx - ay);
    rot[3] = 2 * (xy - az);
    rot[4] = a2 - x2 + y2 - z2;
    rot[5] = 2 * (yz + ax);
    rot[6] = 2 * (zx + ay);
    rot[7] = 2 * (yz - ax);
    rot[8] = a2 - x2 - y2 + z2;
    
    return (1);
}


static void
CenterCoords(double **coords, const int len, const double *weight)
{
    int             i;
    double          xsum, ysum, zsum, wsum;
    double         *x = coords[0], *y = coords[1], *z = coords[2];
    
    xsum = ysum = zsum = 0.0;
    
    if (weight != NULL)
    {
        wsum = 0.0;
        for (i = 0; i < len; ++i)
        {
            xsum += weight[i] * x[i];
            ysum += weight[i] * y[i];
            zsum += weight[i] * z[i];
            
            wsum += weight[i];
        }
        
        xsum /= wsum;
        ysum /= wsum;
        zsum /= wsum;
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            xsum += x[i];
            ysum += y[i];
            zsum += z[i];
        }
        
        xsum /= len;
        ysum /= len;
        zsum /= len;
    }
    
    for (i = 0; i < len; ++i)
    {
        x[i] -= xsum;
        y[i] -= ysum;
        z[i] -= zsum;
    }
}


/* Superposition coords2 onto coords1 -- in other words, coords2 is rotated, coords1 is held fixed */
double
CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight)
{
    double A[9], rmsd;
    /* center the structures -- if precentered you can omit this step */
    CenterCoords(coords1, len, weight);
    CenterCoords(coords2, len, weight);
    
    /* calculate the (weighted) inner product of two structures */
    double E0 = InnerProduct(A, coords1, coords2, len, weight);
    
    /* calculate the RMSD & rotational matrix */
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, len, -1);
    
    return rmsd;
}




/***********************************************************************************************************************************************************************************************************************************************************************************************/
/*
 * Function:  The constructor of EDMD class
 *
 * Parameter: None
 *
 * Return:    None
 */
EDMD::EDMD(void) {
    
    constants.Boltzmann = 0.002;
    constants.timefac = 20.455;
    
    circular = true;
    RNG_Seed = 13579.0;
    
    dt    = 0.002;
    gamma = 0.4;
    tautp = 0.2;
    
    // Convert time-related parameters to internal units
    dt    *= constants.timefac;
    gamma /= constants.timefac;
    tautp *= constants.timefac;
    
    temperature = 300.0;
    
    // Scale factor
    scaled = constants.Boltzmann * temperature;
}


/*
 * Function:  Calculate Ed forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_ED_Forces(Tetrad tetrad, float* ED_Forces, float scaled, int ED_Energy) {

    int i, j;
    double rotmat[9], rmsd, temp;
    double *temp_Forces, *proj;
    
    temp_Forces = new double[3 * tetrad.num_Atoms_In_Tetrad];
    proj = new double[tetrad.num_Evecs];
    
    // Allocate memory for arrays that used in the QCP rotation calculation
    double **avg_Crds  = new double*[3];
    double **crds      = new double*[3];
    double **temp_Crds = new double*[3];
    for (i = 0; i < 3; i++) {
        avg_Crds[i]  = new double[tetrad.num_Atoms_In_Tetrad];
        crds[i]      = new double[tetrad.num_Atoms_In_Tetrad];
        temp_Crds[i] = new double[tetrad.num_Atoms_In_Tetrad];
    }
    
    // Copy data from tetrads
    for (i = 0; i < 3; i++) {
        for (j = 0; j < tetrad.num_Atoms_In_Tetrad; j++) {
            avg_Crds[i][j] = (double) tetrad.avg_Structure[3*i+j];
            crds[i][j]     = (double) tetrad.coordinates[3*i+j];
        }
    }
    
    // Step 1: rotate x into the pcz frame of reference
    rmsd = CalcRMSDRotationalMatrix((double **) avg_Crds, (double **) crds, tetrad.num_Atoms_In_Tetrad, rotmat, NULL);
    for (i = 0; i < tetrad.num_Atoms_In_Tetrad; i++) {
        temp_Crds[0][i] = rotmat[0]*crds[0][i] + rotmat[1]*crds[1][i] + rotmat[2]*crds[2][i];
        temp_Crds[1][i] = rotmat[3]*crds[0][i] + rotmat[4]*crds[1][i] + rotmat[5]*crds[2][i];
        temp_Crds[2][i] = rotmat[6]*crds[0][i] + rotmat[7]*crds[1][i] + rotmat[8]*crds[2][i];
    }
    
    // Step 2: remove average structure, then calculate projections
    for (i = 0; i < tetrad.num_Atoms_In_Tetrad; i++) {
        temp_Crds[0][i] -= avg_Crds[0][i];
        temp_Crds[1][i] -= avg_Crds[1][i];
        temp_Crds[2][i] -= avg_Crds[2][i];
    }
    for(i = 0; i < tetrad.num_Evecs; i++) {
        for(j = 0; j < tetrad.num_Atoms_In_Tetrad; j++) {
            proj[i] += tetrad.eigenvectors[i][j]   * temp_Crds[0][j];
            proj[i] += tetrad.eigenvectors[i][j+1] * temp_Crds[1][j];
            proj[i] += tetrad.eigenvectors[i][j+2] * temp_Crds[2][j];
        }
    }
    
    // Step 3 & Step 4
    for (i = 0; i < tetrad.num_Atoms_In_Tetrad; i++) {
        temp_Crds[0][i] = avg_Crds[0][i];
        temp_Crds[1][i] = avg_Crds[1][i];
        temp_Crds[2][i] = avg_Crds[2][i];
    }
    for (i = 0; i < tetrad.num_Evecs; i++) {
        for (j = 0; j < tetrad.num_Atoms_In_Tetrad; j++) {
            
            // Step 3: re-embed the input coordinates in PC space - a sort of 'shake' procedure.
            //         Ideally this step is not needed, as stuff above should ensure all moves remain in PC subspace...
            temp_Crds[0][j] += tetrad.eigenvectors[i][j]  *proj[i];
            temp_Crds[1][j] += tetrad.eigenvectors[i][j+1]*proj[i];
            temp_Crds[2][j] += tetrad.eigenvectors[i][j+2]*proj[i];
            
            // Step 4: calculate ED forces
            ED_Forces[j]   -= (tetrad.eigenvectors[i][j]  *proj[i]*scaled/tetrad.eigenvalues[i]);
            ED_Forces[j+1] -= (tetrad.eigenvectors[i][j+1]*proj[i]*scaled/tetrad.eigenvalues[i]);
            ED_Forces[j+2] -= (tetrad.eigenvectors[i][j+2]*proj[i]*scaled/tetrad.eigenvalues[i]);
        }
    }
    
    // Step 5: rotate 'shaken' coordinates back into right frame
    for (i = 0; i < tetrad.num_Atoms_In_Tetrad; i++) {
        avg_Crds[0][i] = rotmat[0]*temp_Crds[0][i] + rotmat[1]*temp_Crds[1][i] + rotmat[2]*temp_Crds[2][i];
        avg_Crds[1][i] = rotmat[3]*temp_Crds[0][i] + rotmat[4]*temp_Crds[1][i] + rotmat[5]*temp_Crds[2][i];
        avg_Crds[2][i] = rotmat[6]*temp_Crds[0][i] + rotmat[7]*temp_Crds[1][i] + rotmat[8]*temp_Crds[2][i];
    }
    for (<#initialization#>; <#condition#>; <#increment#>) {
        tetrad.coordinates = avg_Crds;
    }
    
    // Step 6: rotate forces back to original orientation of coordinates
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad;) {
        temp_Forces[i]   = rotmat[0]*ED_Forces[i] + rotmat[1]*ED_Forces[i+1] + rotmat[2]*ED_Forces[i+2];
        temp_Forces[i+1] = rotmat[3]*ED_Forces[i] + rotmat[4]*ED_Forces[i+1] + rotmat[5]*ED_Forces[i+2];
        temp_Forces[i+2] = rotmat[6]*ED_Forces[i] + rotmat[7]*ED_Forces[i+1] + rotmat[8]*ED_Forces[i+2];
        i += 3;
    }
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; i++) {
        ED_Forces[i] = temp_Forces[i];
    }
    
    // Step 7: calculate the 'potential energy' (in units of kT)
    for (i = 0; i < tetrad.num_Evecs; i++) {
        temp += proj[i]*proj[i] / tetrad.eigenvalues[i];
    }
    ED_Forces[ED_Energy] = scaled * 0.5 * temp; // ED Energy
    
    
    for (i = 0; i < 3; i++) {
        delete [] avg_Crds[i];
        delete [] crds[i];
        delete [] temp_Crds[i];
    }
    delete [] avg_Crds;
    delete [] crds;
    delete [] temp_Crds;
    delete [] temp_Forces;
    delete [] proj;
    
    /*
    int i, j, find_Move = 1;
    float temp = 0.0, rmsd = 0.0;
    float *temp_Crds, *temp_Forces, *projections;
    float rotation[3][3], translation[3] = {0.0, 0.0, 0.0};
    
    temp_Crds   = new float[3 * tetrad.num_Atoms_In_Tetrad];
    temp_Forces = new float[3 * tetrad.num_Atoms_In_Tetrad];
    projections = new float[tetrad.num_Evecs];
    
    // Step 1: rotate x into the pcz frame of reference
    //MATFIT(&num_Atoms_In_Tetrad, tetrad.avg_Structure, tetrad.coordinates, rotation, translation, &rmsd, &find_Move);
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; ) {
        temp_Crds[i]   = rotation[0][0]*tetrad.coordinates[i] + rotation[1][0]*tetrad.coordinates[i+1] + rotation[2][0]*tetrad.coordinates[i+2] + translation[0];
        temp_Crds[i+1] = rotation[0][1]*tetrad.coordinates[i] + rotation[1][1]*tetrad.coordinates[i+1] + rotation[2][1]*tetrad.coordinates[i+2] + translation[1];
        temp_Crds[i+2] = rotation[0][2]*tetrad.coordinates[i] + rotation[1][2]*tetrad.coordinates[i+1] + rotation[2][2]*tetrad.coordinates[i+2] + translation[2];
        i = i + 3;
    }
    
    // Step 2: remove average structure, then calculate projections
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; i++) {
        temp_Crds[i] -= tetrad.avg_Structure[i];
    }
    for(i = 0; i < tetrad.num_Evecs; i++) {
        for(j = 0; j < 3 * tetrad.num_Atoms_In_Tetrad; j++) {
            projections[i] += tetrad.eigenvectors[i][j] * temp_Crds[j];
        }
    }
    
    // Step 3 & Step 4
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; i++) {
        temp_Crds[i] = tetrad.avg_Structure[i];
    }
    for (i = 0; i < tetrad.num_Evecs; i++) {
        for (j = 0; j < 3 * tetrad.num_Atoms_In_Tetrad; j++) {
            // Step 3: re-embed the input coordinates in PC space - a sort of 'shake' procedure.
            //         Ideally this step is not needed, as stuff above should ensure all moves remain in PC subspace...
            temp_Crds[j] = temp_Crds[j] + tetrad.eigenvectors[i][j]*projections[i];
            
            // Step 4: calculate ED forces
            ED_Forces[j] = ED_Forces[j] - (tetrad.eigenvectors[i][j]*projections[i]*scaled/tetrad.eigenvalues[i]);
        }
    }
    
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; ) {
        
        // Step 5: rotate 'shaken' coordinates back into right frame
        temp_Crds[i]   -= translation[0];
        temp_Crds[i+1] -= translation[1];
        temp_Crds[i+2] -= translation[2];
        tetrad.coordinates[i]   = rotation[0][0]*temp_Crds[i] + rotation[0][1]*temp_Crds[i+1] + rotation[0][2]*temp_Crds[i+2];
        tetrad.coordinates[i+1] = rotation[1][0]*temp_Crds[i] + rotation[1][1]*temp_Crds[i+1] + rotation[1][2]*temp_Crds[i+2];
        tetrad.coordinates[i+2] = rotation[2][0]*temp_Crds[i] + rotation[2][1]*temp_Crds[i+1] + rotation[2][2]*temp_Crds[i+2];
        
        // Step 6: rotate forces back to original orientation of coordinates
        temp_Forces[i]   = rotation[0][0]*ED_Forces[i] + rotation[0][1]*ED_Forces[i+1] + rotation[0][2]*ED_Forces[i+2];
        temp_Forces[i+1] = rotation[1][0]*ED_Forces[i] + rotation[1][1]*ED_Forces[i+1] + rotation[1][2]*ED_Forces[i+2];
        temp_Forces[i+2] = rotation[2][0]*ED_Forces[i] + rotation[2][1]*ED_Forces[i+1] + rotation[2][2]*ED_Forces[i+2];
        
        i = i + 3;
    }
    for (i = 0; i < 3 * tetrad.num_Atoms_In_Tetrad; i++) {
        ED_Forces[i] = temp_Forces[i];
    }
    
    // Step 7: calculate the 'potential energy' (in units of kT)
    for (i = 0; i < tetrad.num_Evecs; i++) {
        temp += projections[i]*projections[i] / tetrad.eigenvalues[i];
    }
    ED_Forces[ED_Energy] = scaled * 0.5 * temp; // ED Energy
    // *ED_Energy = scaled * 0.5 * temp;
    
    delete []temp_Crds;
    delete []temp_Forces;
    delete []projections;
     */
}


/*
 * Function:  Calculate NB forces
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::calculate_NB_Forces(Tetrad tetrad[], float** NB_Forces, int NB_Energy, int Electrostatic_Energy) {
    
    int i, j;
    float dx, dy, dz, sqdist;
    float a, pair_Force;
    float krep = 100.0;   // krep: soft repulsion constant
    float q;              // q: num_atoms vectors of charges
    float qfac = 332.064; //qfac: electrostatics factor
    float max_Forces = 1.0;
    
    for (i = 0; i < tetrad[0].num_Atoms_In_Tetrad; i++) {
        for (j = 0;  j < tetrad[1].num_Atoms_In_Tetrad; j++) {
            
            dx = tetrad[0].coordinates[3*i]   - tetrad[1].coordinates[3*i];
            dy = tetrad[0].coordinates[3*i+1] - tetrad[1].coordinates[3*i+1];
            dz = tetrad[0].coordinates[3*i+2] - tetrad[1].coordinates[3*i+2];
            sqdist = dx*dx + dy*dy + dz*dz;
            
            // NB energies
            a = max(0.0, 2.0-sqdist);
            q = tetrad[0].abq[3*i+2] * tetrad[1].abq[3*j+2];
            
            NB_Forces[0][NB_Energy] += 0.25 * krep * a * a;
            NB_Forces[1][Electrostatic_Energy] += 0.5 * qfac * q * sqdist;
            // *NB_Energy += 0.25 * krep * a * a;
            // *Electrostatic_Energy += 0.5 * qfac * q * sqdist;
            
            // NB forces
            pair_Force = -2.0 * krep * a - 2.0 * qfac * q / (sqdist * sqdist);
            NB_Forces[0][3*i]   -= dx * pair_Force;
            NB_Forces[0][3*i+1] -= dy * pair_Force;
            NB_Forces[0][3*i+2] -= dz * pair_Force;
            
            NB_Forces[1][3*j]   += dx * pair_Force;
            NB_Forces[1][3*j+1] += dy * pair_Force;
            NB_Forces[1][3*j+2] += dz * pair_Force;
        }
    }
    
    // Clip NB forces
    for (i = 0; i < 3 * tetrad[0].num_Atoms_In_Tetrad; i++) {
        NB_Forces[0][i] = max( max_Forces, NB_Forces[0][i]);
        NB_Forces[0][i] = min(-max_Forces, NB_Forces[0][i]);
    }
    for (i = 0; i < 3 * tetrad[1].num_Atoms_In_Tetrad; i++) {
        NB_Forces[1][i] = max( max_Forces, NB_Forces[1][i]);
        NB_Forces[1][i] = min(-max_Forces, NB_Forces[1][i]);
    }
}


/*
 * Function:  Generate the Gaussian stochastic term. Assuming unitless.
 *
 * Parameter:
 *
 * Return:    None
 */
float EDMD::generate_Stochastic_Term(float tetrad_ID) {
    
    int RNG = RNG_Seed + tetrad_ID;
    
    return RNG;
}


/*
 * Function:  Generate the pair list of tetrads.
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::generate_Pair_Lists(int pair_List[][2], int num_Tetrads, Tetrad* tetrad) {
    
    float mole_Cutoff = 40.0; // Molecular cutoffs
    float atom_Cutoff = 10.0; // Atomic cutoffs
    float mole_Least  = 4.0;  // Molecules less than NB_Cutoff won't have NB ints.
    
    int i, j, k, num_Pairs = 0;
    float r, com[num_Tetrads][3];
    
    // The centre of mass (actually, centre of geom)
    // com(1) = sum(x(1:(3*natoms-2):3))/natoms
    // com(2) = sum(x(2:(3*natoms-1):3))/natoms
    // com(3) = sum(x(3:(3*natoms):3))  /natoms
    for (i = 0; i < num_Tetrads; i++) {
        for (j = 0; j < 3 * tetrad[i].num_Atoms_In_Tetrad; ) {
            com[i][0] += tetrad[i].coordinates[j];
            com[i][1] += tetrad[i].coordinates[j+1];
            com[i][2] += tetrad[i].coordinates[j+2];
            j = j + 3;
        }
        com[i][0] /= tetrad[i].num_Atoms_In_Tetrad;
        com[i][1] /= tetrad[i].num_Atoms_In_Tetrad;
        com[i][2] /= tetrad[i].num_Atoms_In_Tetrad;
    }
    
    // Loop to generate pairlists
    for (i = 0; i < num_Tetrads; i++) {
        for (j = i+1; j < num_Tetrads; j++) {
            r = 0.0;
            
            // If r exceeds mole_Cutoff then no interaction between these two mols
            // r = sum( (com(:,i)-com(:,j)) * (com(:,i)-com(:,j)) )
            for (k = 0; k < 3; k++) {
                r += (com[k][i] - com[k][j]) * (com[k][i] - com[k][j]);
            }
            
            if ((r < mole_Cutoff * mole_Cutoff) && (abs(i-j) > mole_Least) &&
                (abs(i-j) < num_Tetrads - mole_Least)) {
                
                pair_List[num_Pairs][0] = i;
                pair_List[num_Pairs][1] = j;
                
            } else {
                pair_List[num_Pairs][0] = -1;
                pair_List[num_Pairs][1] = -1;
            }
            
            num_Pairs++;
        }
    }
}


/*
 * Function:  Update velocities & Berendsen temperature control
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::update_Velocities(float* velocities, float* ED_Forces, float* NB_Forces, float* masses, int total_Atoms) {
    
    int i;
    float kentical_Energy = 0.0;
    float target_KE;
    float actual_Temperature;
    float tscal;  // Berendsen T-coupling factor
    float gamfac; // Velocity scale factor
    
    gamfac = 1.0 / (1.0 + gamma * dt);
    
    // Simple Langevin dynamics
    for (i = 0; i < 3 * total_Atoms; i++) {
        velocities[i] = (velocities[i] + ED_Forces[i]*dt + NB_Forces[i]*dt/masses[i]) * gamfac;
    }
    
    // Berendsen temperature control
    for (i = 0; i < 3 * total_Atoms; i++) {
        kentical_Energy += 0.5*masses[i]*velocities[i]*velocities[i];
    }
    
    actual_Temperature = kentical_Energy * 2 / (constants.Boltzmann * 3 * total_Atoms);
    
    target_KE = 0.5 * scaled * 3 * total_Atoms;
    
    tscal = sqrt(1.0 + (dt/tautp) * ((target_KE/kentical_Energy) - 1.0));
    
    // Update velocities
    for (i = 0; i < 3 * total_Atoms; i++) {
        velocities[i] = velocities[i] * tscal;
    }
}


/*
 * Function:  Update Coordinates of tetrads
 *
 * Parameter:
 *
 * Return:    None
 */
void EDMD::update_Coordinates(float* coordinates, float* velocities, int total_Atoms) {
    
    for (int i = 0; i < 3 * total_Atoms; i++) {
        coordinates[i] = coordinates[i] + velocities[i]*dt;
    }
}











