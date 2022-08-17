#!/bin/bash

cp $1 $1.tmppoly

file=$1.tmppoly

nbasis=$2

nfields=0

nstates=7

pstates=6

echo -e "#!python
#cython: language_level=3
#cython: extra_compile_args=[\"-Ofast\"]

import numpy as np
import cython

cimport numpy as np

np.import_array()

DTYPE = np.float64

ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
cdef nonaffine_2(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params, np.ndarray[DTYPE_t, ndim=3] psi, np.ndarray[DTYPE_t] mesh):

	cdef int i, j, ele
	cdef int nbox = mesh.shape[0]

	cdef np.ndarray[DTYPE_t] result = np.zeros([26], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t] rhon = np.zeros([nbox], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t] rhop = np.zeros([nbox], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t] rhotot = np.zeros([nbox], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] npot = np.zeros([7,nbox], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t, ndim=2] ppot = np.zeros([6,nbox], dtype=DTYPE)

	cdef np.ndarray[DTYPE_t] angmom = np.zeros([13], dtype=DTYPE)

	cdef double cddr0 = params[6] / 16.00
	cdef double cddr1 = - 1.0/24.0 * params[6] * ( 1.0/2.0 + params[7] )

	cdef double pi = 3.141592653589793

	angmom[0] = 0.5
	angmom[1] = 1.5
	angmom[2] = 0.5
	angmom[3] = 2.5
	angmom[4] = 0.5
	angmom[5] = 1.5
	angmom[6] = 3.5
	angmom[7] = 0.5
	angmom[8] = 0.5
	angmom[9] = 1.5
	angmom[10] = 1.5
	angmom[11] = 2.5
	angmom[12] = 0.5

	rhon[0] = 1e-15
	rhop[0] = 1e-15

	for j in range(7):
		ele = j*2
		for i in range(1,nbox):
			rhon[i] += (2*angmom[j]+1)*(x0[ele]*psi[j,0,i] + x0[ele+1]*psi[j,1,i])**2 \\
			 / (4*pi*mesh[i]**2) + 1e-15

	for j in range(7,13):
		ele = j*2
		for i in range(1,nbox):
			rhop[i] += (2*angmom[j]+1)*(x0[ele]*psi[j,0,i] + x0[ele+1]*psi[j,1,i])**2 \\
			 / (4*pi*mesh[i]**2) + 1e-15

	for j in range(7):
		ele = j*2
		for i in range(1,nbox):
			rhotot[i] = rhon[i] + rhop[i]
			
			npot[j,i] = -((2.0+params[9])*(cddr0-cddr1)*rhotot[i]**(params[9]+1) + \\
			2.0*params[9]*cddr1*(rhon[i]**2 + rhop[i]**2)*rhotot[i]**(params[9]-1) + \\
			4.0*cddr1*rhon[i]*rhotot[i]**params[9])*(x0[ele]*psi[j,0,i] + x0[ele+1]*psi[j,1,i])
	
	for j in range(7,13):
		ele = j*2
		for i in range(1,nbox):
			rhotot[i] = rhon[i] + rhop[i]
			
			ppot[j-7,i] = -((2.0+params[9])*(cddr0-cddr1)*rhotot[i]**(params[9]+1) + \\
			2.0*params[9]*cddr1*(rhon[i]**2 + rhop[i]**2)*rhotot[i]**(params[9]-1) + \\
			4.0*cddr1*rhop[i]*rhotot[i]**params[9] - 1.4399784 * (3/pi)**(1.0/3.0)*rhop[i]**(1.0/3.0)) \\
			* (x0[ele]*psi[j,0,i] + x0[ele+1]*psi[j,1,i])


	result[0] =  np.dot(psi[0,0,:],npot[0,:])
	result[1] =  np.dot(psi[0,1,:],npot[0,:])
	result[2] =  np.dot(psi[1,0,:],npot[1,:])
	result[3] =  np.dot(psi[1,1,:],npot[1,:])
	result[4] =  np.dot(psi[2,0,:],npot[2,:])
	result[5] =  np.dot(psi[2,1,:],npot[2,:])
	result[6] =  np.dot(psi[3,0,:],npot[3,:])
	result[7] =  np.dot(psi[3,1,:],npot[3,:])
	result[8] =  np.dot(psi[4,0,:],npot[4,:])
	result[9] =  np.dot(psi[4,1,:],npot[4,:])
	result[10] = np.dot(psi[5,0,:],npot[5,:])
	result[11] = np.dot(psi[5,1,:],npot[5,:])
	result[12] = np.dot(psi[6,0,:],npot[6,:])
	result[13] = np.dot(psi[6,1,:],npot[6,:])

	result[14] =  np.dot(psi[7,0,:],ppot[0,:])
	result[15] =  np.dot(psi[7,1,:],ppot[0,:])
	result[16] =  np.dot(psi[8,0,:],ppot[1,:])
	result[17] =  np.dot(psi[8,1,:],ppot[1,:])
	result[18] =  np.dot(psi[9,0,:],ppot[2,:])
	result[19] =  np.dot(psi[9,1,:],ppot[2,:])
	result[20] = np.dot(psi[10,0,:],ppot[3,:])
	result[21] = np.dot(psi[10,1,:],ppot[3,:])
	result[22] = np.dot(psi[11,0,:],ppot[4,:])
	result[23] = np.dot(psi[11,1,:],ppot[4,:])
	result[24] = np.dot(psi[12,0,:],ppot[5,:])
	result[25] = np.dot(psi[12,1,:],ppot[5,:])

	return result



@cython.boundscheck(False)
@cython.wraparound(False)
def skyrme_poly_2(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params, np.ndarray[DTYPE_t, ndim=3] psi, np.ndarray[DTYPE_t] mesh):

	cdef np.ndarray[DTYPE_t] result = np.ones([x0.shape[0]], dtype=DTYPE)
	cdef np.ndarray[DTYPE_t] nonaffine = np.zeros([26], dtype=DTYPE)

	nonaffine = nonaffine_2(x0, params, psi, mesh)
" > skyrme_poly.tmppoly
#PsiProtons[1, 1] . FNonAffineProtons[Coefficients, Parameters, 1] +
# {"t0":t0,"x0":x0,"t1":t1,"x1":x1,"t2":t2,"x2":x2,"t3":t3,"x3":x3,"w0":t4p,"sig":1/sigma}
sed -i 's/*\^/e/g' $file 
sed -i 's/\^/**/g' $file 
sed -i 's/0\. +//g' $file
#sed -i 's/ /*/g' $file
sed -i 's/t0/params[0]/g' $file
sed -i 's/x0/params[1]/g' $file
sed -i 's/t1/params[2]/g' $file
sed -i 's/x1/params[3]/g' $file
sed -i 's/t2/params[4]/g' $file
sed -i 's/x2/params[5]/g' $file
sed -i 's/t3/params[6]/g' $file
sed -i 's/x3/params[7]/g' $file
sed -i 's/W/params[8]/g' $file
sed -i 's/sig/params[9]/g' $file

for i in `seq 1 $nstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/NeutronsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file
sed -i "s/PsiNeutrons\[$i, $j\] \. FNonAffineNeutrons\[Coefficients, Parameters, $i\] / nonaffine[$ele] /g" $file

done
done

for i in `seq 1 $pstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + nstates*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/ProtonsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file
sed -i "s/PsiProtons\[$i, $j\] \. FNonAffineProtons\[Coefficients, Parameters, $i\] / nonaffine[$ele] /g" $file

done
done

for i in `seq 1 $nstates`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + pstates*nbasis + nstates*nbasis + (i-1)}')

sed -i "s/EN\[$i\]/x0[$ele]/g" $file

done

for i in `seq 1 $pstates`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + pstates*nbasis + nstates*nbasis + nstates + (i-1)}')

sed -i "s/EP\[$i\]/x0[$ele]/g" $file

done



# sed -i "s/PsiProtons\[[0-9], [0-9]\] \. FNonAffineProtons\[Coefficients, Parameters, [0-9]\] //g" $file
# sed -i "s/PsiNeutrons\[[0-9], [0-9]\] \. FNonAffineNeutrons\[Coefficients, Parameters, [0-9]\] //g" $file

sed -i 's/== 0/ /g' $file
#sed -i 's/{/	result = [/g' $file
sed -i 's/{//g' $file
sed -i 's/}//g' $file
#sed -i 's/\*10\*\*/e/g' $file

totalele=$(awk -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + pstates*nbasis + nstates*nbasis + nstates + pstates}')

sed -i 's/,/\n/g' $file

i=0
while read l; do
  echo -e "\tresult[$i]=$l" >> py.tmppoly
  i=$((i+1))
done < $file

cat skyrme_poly.tmppoly py.tmppoly > skyrme_poly_$2.pyx

echo ""  >> skyrme_poly_$2.pyx
echo "	return result" >> skyrme_poly_$2.pyx

#python3 setup.py build_ext --inplace

#cp skyrme_poly.cpython-39-x86_64-linux-gnu.so ../skyrme_poly.so

rm *.tmppoly