
template < typename TYPE, int DIM >
Evol<TYPE,DIM>* ReadEvol(ifstream &f, int verbosemode=1)
{
    Evol<TYPE,DIM> *Ev = NULL;
    string DefType;
    f >> DefType;
#ifdef LARGEDEF
    if(!DefType.compare("LargeDef"))
        Ev = new LargeDef<TYPE,DIM>(f, verbosemode);
    else
#endif // LARGEDEF
#ifdef LARGEDEF_INITPARAM
    if(!DefType.compare("LargeDef_InitParam"))
        Ev = new LargeDef_InitParam<TYPE,DIM>(f, verbosemode);
    else
#endif // LARGEDEF_INITPARAM
#ifdef LARGEDEFSPEC
        if(!DefType.compare("LargeDefSpec"))
            Ev = new LargeDefSpec<TYPE,DIM>(f, verbosemode);
        else
#endif // LARGEDEFSPEC
#ifdef SMALLDEF
            if(!DefType.compare("SmallDef"))
                Ev = new SmallDef<TYPE,DIM>(f, verbosemode);
            else
#endif // SMALLDEF
#ifdef FREEEVOL
                if(!DefType.compare("FreeEvol"))
                    Ev = new FreeEvol<TYPE,DIM>(f, verbosemode);
                else
#endif // FREEEVOL
				{
                    cout << "Error building Evol from file: unknown type '" << DefType << "'" << endl;
					throw -1;
				}
    return Ev;
}


