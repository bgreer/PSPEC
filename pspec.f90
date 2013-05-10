! PSPEC
! Created May 2012 Ben G.

! This code will read a FITS file, determine the dimensions, allocate enough
!  memory, perform the necessary computations, and finally output a FITS file.


PROGRAM PSPEC
	USE omp_lib
	IMPLICIT NONE
	INCLUDE "fftw3.f"

	CHARACTER(LEN=200) :: infile, outfile, record
	INTEGER :: N, ii, ij, ik, il, blocksize, stat, nfnd
	INTEGER :: offset, nv, g, bitpix, id, nmore, nkeys
	LOGICAL :: anyf, simple, extend
	PARAMETER(N=4)
	INTEGER*8 :: plan
	INTEGER, DIMENSION(3) :: naxes
	DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: inarr, inarr2
	DOUBLE COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: outarr, outarr2
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: power, power_shift
	REAL, DIMENSION(:,:), ALLOCATABLE :: mask2d
	REAL, DIMENSION(:), ALLOCATABLE :: mask1d
	REAL :: apod_min, apod_min_t, r, cx, cy, x, y, x0, y0, dtheta, tht, pi
	REAL :: mapscale, dk, dnu
	REAL, DIMENSION(:), ALLOCATABLE :: buff
	INTEGER :: ntheta, nk, ntheta_sub, npix, pixloc, reducefactor
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: power_cart, power_sub
	DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: linein, lineout, subin, subout
	DOUBLE PRECISION :: norm

	pi = ACOS(-1D0)

	! Check number of command line arguments
	IF (IARGC() .lt. 1) THEN
		PRINT*, "Usage: ./pspec infile outfile"
		STOP
	ENDIF

	! Read command line arguments
	CALL getarg(1,infile)
	CALL getarg(2,outfile)

	!$OMP PARALLEL
	id = OMP_GET_THREAD_NUM()
	IF (id .eq. 0) THEN
		print*, "Number of threads: ", OMP_GET_NUM_THREADS()
	ENDIF
	!$OMP END PARALLEL

	stat = 0
	! Open input FITS file
	CALL FTOPEN(30,infile,0,blocksize,stat)
	! Get dimensions of data cube
	CALL FTGKNJ(30,'NAXIS',1,3,naxes,nfnd,stat)
	PRINT*, "Input data cube has dimensions: ",naxes
	npix = naxes(1)

	ALLOCATE(inarr(naxes(1),naxes(2),naxes(3)))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READ INPUT

	! Read input FITS into inarr
	PRINT*, "Reading input file into memory.."
	ALLOCATE(buff(naxes(1)))
	nv = -999
	g = 1
	offset = 1
	buff = 0
!	CALL FTG3DE(30,g,nv,naxes(3),naxes(2),naxes(3),naxes(2),naxes(1),inarr,anyf,stat)
	DO ii=1,naxes(3)
		DO ij=1,naxes(2)
			CALL FTGPVE(30,g,offset,naxes(1),nv,buff,anyf,stat)
			inarr(:,ij,ii) = buff
			offset = offset + naxes(1)
		ENDDO
	ENDDO
	DEALLOCATE(buff)
	PRINT*, " Done."

!	PRINT*, inarr(100,100,1), 151.758

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BACKSUBTR

	!$OMP PARALLEL DO PRIVATE(ii)
	DO ii=1,naxes(3)
		inarr(:,:,ii) = inarr(:,:,ii) - SUM(inarr(:,:,ii)) / &
			(float(naxes(1))*float(naxes(2)))
	ENDDO
	!$OMP END PARALLEL DO
	PRINT*, " Done."

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! APODIZE

	! Create Apodization Masks
	PRINT*, "Apodizing.."
	ALLOCATE(mask2d(naxes(1),naxes(2)))
	ALLOCATE(mask1d(naxes(3)))
	mask1d = 1.0
	mask2d = 1.0
	apod_min = 0.9375
	apod_min_t = 0.96875
	cx = naxes(1)/2.+0.5
	cy = naxes(2)/2.+0.5
	DO ii=1,naxes(3) ! Time
		! missing a factor of 2?
		IF (FLOAT(ii-1)/FLOAT(naxes(3)) .LT. 1.0-apod_min_t) THEN
			mask1d(ii) = (1.-(1.-(FLOAT(ii-1)/FLOAT(naxes(3)))/(1.0-apod_min_t))**2.)**2.
		ENDIF
		IF (FLOAT(ii-1)/FLOAT(naxes(3)) .GT. apod_min_t) THEN
			mask1d(ii) = (1.-(1.-(FLOAT(naxes(3)-ii)/FLOAT(naxes(3)))/(1.-apod_min_t))**2.)**2.
		ENDIF
	ENDDO
	DO ii=1,naxes(2) ! Space
		DO ij=1,naxes(1)
			r = MIN(SQRT((ij-cx)**2.+(ii-cy)**2.) / cx, 1.0)
			IF (r .GT. apod_min) THEN
				mask2d(ii,ij) = (1.-((r-apod_min)/(1.-apod_min))**2.)**2.
			ENDIF
		ENDDO
	ENDDO

	! Apply Apodization to Data
	!$OMP PARALLEL DO PRIVATE(ii)
	DO ii=1,naxes(3)
		inarr(:,:,ii) = inarr(:,:,ii)*mask1d(ii)*mask2d
	ENDDO
	!$OMP END PARALLEL DO
	PRINT*, " Done."
!	PRINT*, inarr(100,100,10), -3.184

	! Free masks
	DEALLOCATE(mask2d)
	DEALLOCATE(mask1d)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOURIER TRANSFORM

	! transform because i'm lazy
	PRINT*, "Transforming array.."
	ALLOCATE(inarr2(naxes(3), naxes(2), naxes(1)))
	! copy inarr to inarr2
	!$OMP PARALLEL DO PRIVATE(ii,ij)
	DO ii=1,naxes(3)
		DO ij=1,naxes(1)
			inarr2(ii,:,ij) = inarr(ij,:,ii)
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO

	DEALLOCATE(inarr)
	ALLOCATE(outarr2(naxes(3)/2+1, naxes(2), naxes(1)))
	PRINT*, " Done."
!	PRINT*, inarr2(10,100,100), -3.184

	plan = 0
	! Set up FFT of inarr->outarr
	CALL DFFTW_PLAN_DFT_R2C_3D(plan,naxes(3),naxes(2),naxes(1),inarr2,outarr2,FFTW_ESTIMATE)
	PRINT*, "Performing Fourier Transform.."
	CALL DFFTW_EXECUTE(plan)
	PRINT*, " Done."
	CALL DFFTW_DESTROY_PLAN(plan)

	DEALLOCATE(inarr2)
!	PRINT*, outarr2(10,100,100), 13499424.5, 2400790.3

	PRINT*, "Transforming array again.."
	! copy outarr2 to outarr

	ALLOCATE(outarr(naxes(1),naxes(2),naxes(3)/2+1))
	!$OMP PARALLEL DO PRIVATE(ii,ij)
	DO ii=1,naxes(3)/2+1
		DO ij=1,naxes(1)
			outarr(ij,:,ii) = outarr2(ii,:,ij)
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO
	DEALLOCATE(outarr2)
	PRINT*, " Done."

!	PRINT*, outarr(100,100,10), 13499424.5, 2400790.3

	ALLOCATE(power(naxes(1),naxes(2),naxes(3)/2+1))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! POWER SPECTRUM

	! Process result
	PRINT*, "Computing Power Spectrum.."
	!$OMP PARALLEL DO PRIVATE(ii)
	DO ii=1,naxes(3)/2+1
		power(:,:,ii) = 0.5*LOG(REAL(outarr(:,:,ii))**2. + &
			AIMAG(outarr(:,:,ii))**2.)
	ENDDO
	!$OMP END PARALLEL DO
!	DO ii=1,naxes(3)/2+1
!		DO ij=1,naxes(1)/2+1
!			power(ij,:,ii) = 0.5*LOG(REAL(outarr(ij,:,ii))**2. + &
!				AIMAG(outarr(ij,:,ii))**2.)
!		ENDDO
!		DO ij=1,naxes(1)/2-1
!			power(naxes(1)/2+ij+1,:,ii) = 0.5*LOG(REAL(outarr(naxes(1)/2-ij+2,:,naxes(3)-ii+1))**2. + &
!				AIMAG(outarr(naxes(1)/2-ij+2,:,naxes(3)-ii+1))**2.)
!		ENDDO
!	ENDDO
	PRINT*, " Done."
	DEALLOCATE(outarr)

!	PRINT*, power(100,100,10), 16.43373

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CART TO POLAR

	! Shift pixels
	PRINT*, "Shifting Pixels.."
	ALLOCATE(power_shift(naxes(1),naxes(2),naxes(3)/2+1))
	power_shift(naxes(1)/2+1:naxes(1),:,:) = power(1:naxes(1)/2,:,:)
	power_shift(1:naxes(1)/2,:,:) = power(naxes(1)/2+1:naxes(1),:,:)
	power(:,naxes(2)/2+1:naxes(2),:) = power_shift(:,1:naxes(2)/2,:)
	power(:,1:naxes(2)/2,:) = power_shift(:,naxes(2)/2+1:naxes(2),:)
	DEALLOCATE(power_shift)
	PRINT*, " Done."

	PRINT*, "Converting to Polar Coordinates.."
	nk = naxes(1)/2
	ntheta = 200
	SELECT CASE (naxes(1))
		CASE (768) ! 32 degrees
			ntheta = 512
		CASE (384) ! 16 degree
			ntheta = 256
		CASE (192) ! 8 degree
			ntheta = 128
		CASE (96) ! 4 degree
			ntheta = 64
		CASE (48) ! 2 degree
			ntheta = 32
	END SELECT
	ALLOCATE(power_cart(ntheta,nk,naxes(3)/2))
	x0 = FLOAT(naxes(1)/2)+1.
	y0 = FLOAT(naxes(2)/2)+1.
	dtheta = 2.0*pi/FLOAT(ntheta)
	!$OMP PARALLEL DO PRIVATE(ii,ij,ik,tht,x,y)
	DO ii=1,naxes(3)/2
		DO ij=1,nk
			DO ik=1,ntheta
				tht = (ik-1)*dtheta
				x = FLOAT(ij)*COS(tht)+x0
				y = FLOAT(ij)*SIN(tht)+y0
				power_cart(ik,ij,ii) = EXP(INTERP(x,y,ii))
			ENDDO
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO
	DEALLOCATE(power)
	PRINT*, " Done."

!	PRINT*, power_cart(1, 20, 10), 1.2367e7

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOURIER SUBSAMPLE

!	ntheta_sub = ntheta
	ntheta_sub = INT(MAX(ntheta/8,16))
	reducefactor = INT(ntheta/ntheta_sub)
	PRINT*, "Subsampling from ",ntheta," to ",ntheta_sub
	ALLOCATE(power_sub(ntheta_sub,nk,naxes(3)/2))
	ALLOCATE(linein(ntheta))
	ALLOCATE(lineout(ntheta))
	ALLOCATE(subin(ntheta))
	ALLOCATE(subout(ntheta))
	DO ii=1,naxes(3)/2
		DO ij=1,nk
			! copy spectrum -> linein
			linein = CMPLX(power_cart(:,ij,ii),0.0)
			! ft linein -> lineout
			CALL DFFTW_PLAN_DFT_1D(plan,ntheta,linein,lineout,FFTW_FORWARD,FFTW_ESTIMATE)
			CALL DFFTW_EXECUTE_DFT(plan, linein, lineout)
			CALL DFFTW_DESTROY_PLAN(plan)
			! copy lineout -> subin
!			subin = lineout(1:ntheta_sub)/FLOAT(ntheta)
			subin = lineout
			subin(ntheta_sub+1:ntheta-ntheta_sub+1) = CMPLX(0.0,0.0)
			! ft subin -> subout
			CALL DFFTW_PLAN_DFT_1D(plan,ntheta,subin,subout,FFTW_BACKWARD,FFTW_ESTIMATE)
			CALL DFFTW_EXECUTE_DFT(plan, subin, subout)
            CALL DFFTW_DESTROY_PLAN(plan)
			! copy subout > spectrum
!			power_sub(:,ij,ii) = SQRT(REAL(subout)**2.+AIMAG(subout)**2.)
!			power_sub(:,ij,ii) = REAL(subout)
			DO ik=1,ntheta_sub
				power_sub(ik,ij,ii) = REAL(subout((ik-1)*reducefactor+1))
			ENDDO
			!power_sub(:,ij,ii) = power_cart(:,ij,ii)
		ENDDO
	ENDDO
	DEALLOCATE(linein)
	DEALLOCATE(lineout)
	DEALLOCATE(subin)
	DEALLOCATE(subout)
	DEALLOCATE(power_cart)
	PRINT*, " Done."

!	PRINT*, power_sub(1, 20, 10), 4.2708e9

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT

	! Output to FITS file
	PRINT*, "Writing Results to File.."
	dnu = 1e6/(45.*naxes(3))
	CALL FTINIT(20,outfile,blocksize,stat)
	simple=.TRUE.
	bitpix=-32
	extend=.TRUE.
	naxes(3) = naxes(3) / 2
	naxes(2) = nk
	naxes(1) = ntheta_sub
	CALL FTPHPR(20,simple,bitpix,3,naxes,0,1,extend,stat)

	! Copy old header info
	CALL FTRDEF(20,stat)
	CALL FTGHSP(30,nkeys,nmore,stat)
	DO ii=7,nkeys
	   CALL FTGREC(30,ii,record,stat)
	   CALL FTPREC(20,record,stat)
	ENDDO
	CALL FTGKYE(30,"MAPSCALE",mapscale,record,stat)
	dk = 360./(mapscale*npix*696.0)
	PRINT*, " Adding header keys:"
	PRINT*, "   DELTA_NU=",dnu
	PRINT*, "   DELTA_K=",dk
	CALL FTPKYE(20,"DELTA_NU",dnu,6,"",stat)
	CALL FTPKYE(20,"DELTA_K",dk,6,"",stat)
	CALL FTCLOS(30,stat)
	g = 1
	offset = 1
!	CALL FTP3DE(20,g,naxes(1),naxes(2),naxes(1),naxes(2),naxes(3),power_sub,stat)
	DO ii=1,naxes(3)
		DO ij=1,naxes(2)
			CALL FTPPRE(20,g,offset,naxes(1),power_sub(:,ij,ii),stat)
			offset = offset + naxes(1)
		ENDDO
	ENDDO
	CALL FTCLOS(20,stat)
	PRINT*, " Done."

	! Deallocate output
	DEALLOCATE(power_sub)



	CONTAINS

	! Subroutine to interpolate
	REAL FUNCTION INTERP(x,y,findex)
		IMPLICIT NONE
		INTEGER :: ii, ij, findex
		REAL :: x, y, weight, temp, res
		! Coordinates of surrounding points
		INTEGER, DIMENSION(4) :: ix, iy
		! Imaginary data points at y=const
		REAL, DIMENSION(4) :: r

		IF (x .LT. 1.0 .OR. x .GT. naxes(1)-1) THEN
			INTERP = 0.0
			RETURN
		ENDIF
		
		IF (y .LT. 1.0 .OR. y .GT. naxes(2)-1) THEN
			INTERP = 0.0
			RETURN
		ENDIF

		! For bicubic convolution interpolation, need nearest 2 points
		!  in each direction. Find them in O(1) time.
		ix(2) = MAX(FLOOR(x), 1)
		iy(2) = MAX(FLOOR(y), 1)
		ix(3) = MIN(ix(2)+1, naxes(1))
		iy(3) = MIN(iy(2)+1, naxes(2))
		ix(1) = MAX(ix(2)-1,1)
		iy(1) = MAX(iy(2)-1,1)
		ix(4) = MIN(ix(2)+2,naxes(1))
		iy(4) = MIN(iy(2)+2,naxes(2))

		! Interpolate at current y
		r = 0.0
		DO ii=1,4 ! x-index
			weight = 0.0
			DO ij=1,4 ! sum over y
				temp = KERNEL(iy(ij)-y)
				weight = weight + temp
				r(ii) = r(ii) + temp*power(ix(ii),iy(ij),findex)
			ENDDO
			r(ii) = r(ii) / weight
		ENDDO

		! Interpolate at current x
		res = 0.0
		weight = 0.0
		DO ii=1,4
			temp = KERNEL(ix(ii)-x)
			weight = weight + temp
			res = res + temp*r(ii)
		ENDDO

		res = res / weight
		INTERP = res
	END FUNCTION INTERP

	REAL FUNCTION KERNEL (x)
		IMPLICIT NONE
		REAL :: x, res
		x = ABS(x)
		IF (x .le. 1.0) THEN
			res = 1.0 + x*x*(-2.5 + 1.5*x)
		ELSE IF (x .lt. 2.0) THEN
			res = 2.0 + x*(-4.0 + x*(2.5 - 0.5*x))
		ELSE
			res = 0.0
		ENDIF
		KERNEL = res
	END FUNCTION KERNEL
END PROGRAM
