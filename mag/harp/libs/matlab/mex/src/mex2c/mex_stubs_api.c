/*
 * mex_stubs_api.c
 * 
 * These are stubs for mex functions that use *only* the mx* API.
 * They are used in both mex_stubs.c and pymex_stubs.c, 
 * which #include this file.
 *
 * Since these use only the API, not the mxArray structure's internals, 
 * they can be shared between files.
 *
 * We tried to move the functions in groups, e.g. the mxIsXXX
 * functions are all pretty trivial.
 *
 * Michael Turmon, 2008
 */

/* no #includes needed */

/* remember, the internals of the mxArray structure will be different
 * here depending on compiler options in the parent file.
 */

const char *mxGetClassName(const mxArray *pa) {
  switch (mxGetClassID(pa)) {
  case mxSPARSE_CLASS:
    return("sparse");  /* obsolete in r13 and up */
  case mxDOUBLE_CLASS:
    return("double"); 
  case mxSINGLE_CLASS:
    return("single"); 
  case mxINT8_CLASS:
    return("int8"); 
  case mxUINT8_CLASS:
    return("uint8"); 
  case mxINT16_CLASS:
    return("int16"); 
  case mxUINT16_CLASS:
    return("uint16"); 
  case mxINT32_CLASS:
    return("int32"); 
  case mxUINT32_CLASS:
    return("uint32"); 
  case mxINT64_CLASS:
    return("int64"); 
  case mxUINT64_CLASS:
    return("uint64"); 
  case mxCHAR_CLASS: 
    return("char"); 
  case mxLOGICAL_CLASS:
    return("logical"); 
  case mxFUNCTION_CLASS:
    return("function handle");
  case mxOBJECT_CLASS:
    return("object");
  case mxOPAQUE_CLASS:
    return("opaque"); 
  case mxUNKNOWN_CLASS:
    return("unknown");
  default:
    return("unrecognized");
  }
}


size_t
mxGetElementSize(const mxArray *pa) {
  switch(mxGetClassID(pa)) {
  case mxSTRUCT_CLASS:
  case mxOBJECT_CLASS:
    /* believe I need to multiply by number of fields */
    return sizeof(mxArray *); /* the pointer stores them, see doc */
  case mxCELL_CLASS:
    return sizeof(mxArray *); /* the pointer stores them, see doc */
  case mxSPARSE_CLASS:
    return sizeof(double);    /* note: this is obsolete in r13 and up */
  case mxDOUBLE_CLASS:
    return sizeof(double);
  case mxSINGLE_CLASS:
    return sizeof(float);
  case mxCHAR_CLASS: 
    return sizeof(mxChar);    /* mxChar's should be 16-bit, but be sure */
  case mxLOGICAL_CLASS: 
    return sizeof(mxLogical);
  case mxINT8_CLASS:
  case mxUINT8_CLASS:
    return 1;
  case mxINT16_CLASS:
  case mxUINT16_CLASS:
    return 2;
  case mxINT32_CLASS:
  case mxUINT32_CLASS:
    return 4;
  case mxINT64_CLASS:
  case mxUINT64_CLASS:
    return 8;
  case mxFUNCTION_CLASS:
  case mxOPAQUE_CLASS:
  default:
    return 0; /* should perhaps raise an error */
  }
}

bool
mxIsFromGlobalWS(const mxArray *pa)
{
  return(0); /* does not have an interpretation now */
}

bool
mxIsEmpty(const mxArray *pm)
{
  return(mxGetNumberOfElements(pm) == 0);
}

bool
mxIsNumeric(const mxArray *pm)
{
  switch (mxGetClassID(pm)) {
  case mxSPARSE_CLASS:  /* obsolete in r13 and up */
  case mxDOUBLE_CLASS:
  case mxSINGLE_CLASS:
  case mxINT8_CLASS:
  case mxUINT8_CLASS:
  case mxINT16_CLASS:
  case mxUINT16_CLASS:
  case mxINT32_CLASS:
  case mxUINT32_CLASS:
  case mxINT64_CLASS:
  case mxUINT64_CLASS:
    return(1);
    /*NOTREACHED*/
    break;
  default:
    /* includes logical class */
    return(0);
    /*NOTREACHED*/
    break;
  }
}

bool
mxIsDouble(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxDOUBLE_CLASS);
}

bool
mxIsChar(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxCHAR_CLASS); 
}

bool
mxIsCell(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxCELL_CLASS); 
}

bool
mxIsStruct(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxSTRUCT_CLASS); 
}

bool
mxIsFunctionHandle(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxFUNCTION_CLASS); 
}

bool
mxIsObject(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxOBJECT_CLASS); 
}

bool
mxIsOpaque(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxOPAQUE_CLASS); 
}

bool
mxIsClass(const mxArray *pm, const char *name)
{
  mex2c_unimplemented("mxIsClass");
  return (0); /* not supported */
}

bool
mxIsSingle(const mxArray *pm)
{
  return (mxGetClassID(pm) == mxSINGLE_CLASS); 
}

bool
mxIsLogical(const mxArray *pa)
{
  /* the semantics of "logical" was changed in r13; this library follows
   * the new semantics where logical is a separate class like double or char */
  return (mxGetClassID(pa) == mxLOGICAL_CLASS); 
}

bool
mxIsLogicalScalar(const mxArray *pa)
{
  return mxIsLogical(pa) && (mxGetNumberOfElements(pa) == 1);
}

bool
mxIsLogicalScalarTrue(const mxArray *pa)
{
  return mxIsLogicalScalar(pa) && mxGetLogicals(pa)[0];
}

bool
mxIsInt8(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxINT8_CLASS); 
}

bool
mxIsUint8(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxUINT8_CLASS); 
}

bool
mxIsInt16(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxINT16_CLASS); 
}

bool
mxIsUint16(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxUINT16_CLASS); 

}

bool
mxIsInt32(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxINT32_CLASS); 
}

bool
mxIsUint32(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxUINT32_CLASS); 
}

bool
mxIsInt64(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxINT64_CLASS); 
}

bool
mxIsUint64(const mxArray *pa)
{
  return (mxGetClassID(pa) == mxUINT64_CLASS); 
}



