#ifndef TRIANGULATION
#define TRIANGULATION

template < typename TYPE >
class Triangulation
{
        typedef TinyVector<TYPE,3> Vect;
        typedef Array<Vect,1> ArrVect;
        typedef TinyVector<int,3> Face;
        typedef Array<Face,1> ArrFace;

    public:

        ArrVect Vertices;
        ArrFace Faces;

        Triangulation() { }

        Triangulation(int n)
        {
            Init(n);
        }

        void Init(int n)
        {
            Faces.resize(Range(1,n));
            Faces.reindexSelf(1);
        }

        Triangulation(ArrVect &points, Range R = Range::all())
        {
            Reference(points,R);
        }

        Triangulation(Triangulation &T) : Vertices(Range(1,T.Vertices.rows())), Faces(Range(1,T.Faces.rows()))
        {
            Vertices = T.Vertices.copy();
            Vertices.reindexSelf(1);
            Faces = T.Faces.copy();
            Faces.reindexSelf(1);
        }

        void operator=(Triangulation &dm)
        {
            Init(dm.Vertices.rows());
            Vertices = dm.Vertices.copy();
            Vertices.reindexSelf(1);
            Faces = dm.Faces.copy();
            Faces.reindexSelf(1);
        }

        void Reference(ArrVect &points, Range R = Range::all())
        {
            Vertices.reference(points(R));
            Vertices.reindexSelf(1);
        }

        int NumVertices()
        {
            return Vertices.rows();
        }

        int NumFaces()
        {
            return Faces.rows();
        }

        ArrVect ComputeCenters()
        {
            ArrVect Centers(Range(1,Faces.rows()));
            TYPE oot = (TYPE)1.0/(TYPE)3.0;
            for(int i=1; i<F.rows()+1; i++)
                Centers(i) = oot*(Vertices(Faces(i)(0))+Vertices(Faces(i)(1))+Vertices(Faces(i)(2)));
            return Centers;
        }

        ArrVect ComputeNormals()
        {
            ArrVect Normals(Range(1,Faces.rows()));
            Vect u,v;
            for(int i=1; i<Faces.rows()+1; i++)
            {
                u = Vertices(Faces(i)(1))-Vertices(Faces(i)(0));
                v = Vertices(Faces(i)(2))-Vertices(Faces(i)(0));
                u = cross(u,v);
                Normals(i) = 0.5 * u;
            }
            return Normals;
        }

        virtual ~Triangulation() {};
};

#endif // TRIANGULATION
