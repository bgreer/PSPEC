
! Contains global variables, tracking settings, etc.

MODULE ParseInput

	IMPLICIT NONE

	INTEGER :: verbose
	LOGICAL :: dofilter
	CHARACTER(LEN=200) :: infile, outfile

CONTAINS

	SUBROUTINE Usage()
		IMPLICIT NONE
		PRINT*, "Minimal Usage: ./frack"
		PRINT*, "Other Options:"
		PRINT*, " -v      Verbose Mode"
		PRINT*, " -vv     Very Verbose Mode"
		PRINT*, " -f      Perform Fourier Filter"
	END SUBROUTINE Usage

	! Note the lack of IMPLICIT NONE
	! Poor coding? Maybe.
	SUBROUTINE ReadCommandLine()
		INTEGER :: ii, argcount, filenumber
		REAL :: tempreal
		CHARACTER(LEN=200) :: strbuffer

		argcount = IARGC()
		filenumber = 1

		! Read command line arguments
		DO ii=1,argcount
			CALL getarg(ii,strbuffer)

			IF (strbuffer .EQ. "--help" .OR. strbuffer .EQ. "-h") THEN
				CALL Usage()
				STOP
			ENDIF

			! verbose mode
			IF (strbuffer .EQ. "-v") THEN
				verbose = 1

			ELSEIF (strbuffer .EQ. "-vv") THEN
				verbose = 2

			ELSEIF (strbuffer .EQ. "-f") THEN
				dofilter = .TRUE.

			ELSE
				! input or output file name
				IF (filenumber == 1) THEN
					CALL getarg(ii,infile)
					filenumber = 2
				ELSEIF (filenumber == 2) THEN
					CALL getarg(ii,outfile)
					filenumber = 3
				ELSE
					PRINT*, "Invalid command-line argument: ", strbuffer
				ENDIF
			ENDIF
		ENDDO
	END SUBROUTINE ReadCommandLine



	! Set the default values of various run parameters
	! It is best to call this before reading in any conflicting information..
	SUBROUTINE SetDefaults()
		verbose = 0
		dofilter = .FALSE.
	END SUBROUTINE SetDefaults


END MODULE
