C *************************************************************************
      BLOCK DATA

C Module name
        CHARACTER*127 MNAME
        COMMON / MNAME / MNAME
        DATA MNAME / 'ringfit_ssw' /

        CHARACTER*927 MARGS(128)
        COMMON / MARGS / MARGS
	DATA MARGS (1) / 'string,record,
     +unspecified,
     +Input data record'/ 
        DATA MARGS (2) / 'flag,f,unspecified,
     +Run fourier filter'/ 
        DATA MARGS (3) / 'end' /
      END BLOCK DATA

      INTEGER FUNCTION DOIT()
	USE FDRMS 
	integer params, status
        character*256 inds
        INTEGER filtr
        INTEGER nrecs
        INTEGER naxis
        CHARACTER*256 umwelt, spectra, rec, pspec, parmseg
        CHARACTER*256 arrpowxy
        CHARACTER, ALLOCATABLE, TARGET :: datav(:)
	CHARACTER*256 pathname

        CHARACTER*256 segstr, oseries, seqnostr, filo
        DOUBLE PRECISION, ALLOCATABLE, TARGET :: powxy(:)
        INTEGER axis(DRMS_MAXRANK)

	CALL CPGETHANDLE (params)
	CALL CPGETSTR (params, 'record', inds, status)
	CALL CPFLAGSET (params, 'f', filtr, status)

	umwelt = F_DRMS_ENV_HANDLE()

        spectra = F_DRMS_OPEN_RECORDS(umwelt, inds, status)
        nrecs = F_DRMS_RECORDSET_GETNRECS(spectra)
        if (nrecs.NE.1) THEN
           status = F_DRMS_CLOSE_RECORDS(spectra, DRMS_FREE_RECORD)
           DOIT = 0
        ENDIF

        rec = F_DRMS_RECORDSET_GETREC(spectra, 0)
        segstr = 'power'
        pspec = F_DRMS_SEGMENT_LOOKUP(rec, segstr)
        
C       Should check for spec null, but can't
        naxis = F_DRMS_SEGMENT_GETNAXIS(pspec)
        if (naxis.NE.3) THEN
           status = F_DRMS_CLOSE_RECORDS(spectra, DRMS_FREE_RECORD)
           DOIT = 0
        ENDIF

        arrpowxy = F_DRMS_SEGMENT_READ(pspec, DRMS_TYPE_DOUBLE, status)
        CALL F_GET_ARRAY_DATA_DOUBLE(powxy, arrpowxy)
        CALL F_GET_ARRAY_AXIS (axis, arrpowxy)

C       CALL ringanalysis (poxwy, axis(1), axis(2), axis(3))
        oseries = 'hmi.ssw_ringfits16'
        rec = F_DRMS_CREATE_RECORD(umwelt, oseries,
     2                                 DRMS_PERMANENT, status)

        seqnostr = 'seqno'
        status = F_DRMS_SETKEY_INT(rec, seqnostr, 0)

        segstr = 'fit'
        parmseg = F_DRMS_SEGMENT_LOOKUP(rec, segstr)
        filo = 'filename'
        status = F_DRMS_SEGMENT_WRITE_FROM_FILE(parmseg, filo)
        status = F_DRMS_CLOSE_RECORD(rec, DRMS_INSERT_RECORD)
        status = F_DRMS_CLOSE_RECORDS(spectra, DRMS_FREE_RECORD)

	DOIT = 0
	RETURN
      end
