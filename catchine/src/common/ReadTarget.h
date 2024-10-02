#ifdef TARGET


template < typename TYPE, int DIM >
Array<Target<TYPE,DIM>*,1>* ReadArrTargets(ifstream &f, int verbosemode=1)
{
    int Ntargets;
    f >> Ntargets;
    // Create the array of targets
    Array<Target<TYPE,DIM>*,1>* Tgt = new Array<Target<TYPE,DIM>*,1>(Range(1,Ntargets));
    string TargetType;
    for(int i=1; i<Ntargets+1; i++)
    {
        f >> TargetType;
#ifdef MEASURE
        if(!TargetType.compare("Measure"))
            (*Tgt)(i) = new Measure<TYPE,DIM>(f, verbosemode);
        else
#endif // MEASURE
#ifdef LANDMARKS
            if(!TargetType.compare("Landmarks"))
                (*Tgt)(i) = new Landmarks<TYPE,DIM>(f, verbosemode);
            else
#endif // LANDMARKS
#ifdef LANDMARKS_NORM1
            if(!TargetType.compare("Landmarks_Norm1"))
                (*Tgt)(i) = new Landmarks_Norm1<TYPE,DIM>(f, verbosemode);
            else
#endif // LANDMARKS_NORM1
#ifdef SURFCURR
                if(!TargetType.compare("SurfCurr"))
                    (*Tgt)(i) = new SurfCurr<TYPE>(f, verbosemode);
                else
#endif // SURFCURR  
#ifdef CURVECURR
                    if(!TargetType.compare("CurveCurr"))
                        (*Tgt)(i) = new CurveCurr<TYPE,DIM>(f, verbosemode);
                    else
#endif // CURVECURR  
#ifdef CURVECYCLE
                        if(!TargetType.compare("CurveCycle"))
                            (*Tgt)(i) = new CurveCycleTarget<TYPE,DIM>(f, verbosemode);
                        else
#endif // CURVECYCLE  
#ifdef CURVEACC
                            if(!TargetType.compare("CurveAcc"))
                                (*Tgt)(i) = new CurveAcc<TYPE>(f, verbosemode);
                            else
#endif // CURVEACC
#ifdef L2IMAGE
                                if(!TargetType.compare("L2Image"))
                                    (*Tgt)(i) = new L2Image<TYPE>(f, verbosemode);
                                else
#endif // L2IMAGE
								{
                                    cout << "Error reading Target: unknown type '" << TargetType << "'" << endl;
									throw -1;
								}
    }
    return Tgt;
}


#endif // TARGET

