#ifdef OPTIM

template < typename TYPE >
Optim<TYPE>* ReadOptim(ifstream &f, int verbosemode)
{
    Optim<TYPE> *Opt = NULL;
    string Tag;
    f >> Tag;
#ifdef FIXEDGRADDESCENT
    if(!Tag.compare("FixedGradDescent"))
        Opt = new FixedGradDescent<TYPE>(f, verbosemode);
    else
#endif // FIXEDGRADDESCENT
#ifdef ADAPTGRADDESCENT
		if(!Tag.compare("AdaptGradDescent"))
			Opt = new AdaptGradDescent<TYPE>(f, verbosemode);
		else
#endif // ADAPTGRADDESCENT
#ifdef LBFGS_LIBLBFGS
			if(!Tag.compare("Lbfgs_liblbfgs"))
				Opt = new Lbfgs_liblbfgs<TYPE>(f, verbosemode);
			else
#endif // LBFGS_LIBLBFGS
#ifdef LBFGS_LIBLBFGS_QUENTIN
			if(!Tag.compare("Lbfgs_liblbfgs_Quentin"))
				Opt = new Lbfgs_liblbfgs_Quentin<TYPE>(f, verbosemode);
			else
#endif // LBFGS_LIBLBFGS_QUENTIN
#ifdef OPTIM_DLIB
			if(!Tag.compare("Optim_Dlib"))
				Opt = new Optim_Dlib<TYPE>(f, verbosemode);
			else
#endif // OPTIM_DLIB
			{
        cout << "Error constructing Optim from file: unknown type '" << Tag << "'." << endl;
	throw -1;
    }
    return Opt;
}

#endif // OPTIM
