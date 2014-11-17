MyProject
=========

To develop a computer system forming the total sample of the general population subjects when analyzing the results of mass testing.














Код построения выборочных совокупностей: случайной, механической, типической, серийной и двухступенчатой.
namespace SampleStatistics
{
    public class DataStat
    {
        public static void TypicalSampleP(int N, int[,] r, int n, int[] y)
        {
            int temp0, temp1;
            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            List<int> Ni = new List<int>();
            int k = 1;
            for (int i = 0; i < N - 1; i++)
            {
                if (r[i, 0] == r[i + 1, 0]) k++;
                else
                {
                    Ni.Add(k); k = 1;
                }
            }


            int ni, start = 0,start2=0;
            Random rd = new Random();
            for (int j = 0; j < Ni.Count(); j++)
            {
                ni = (int)Math.Ceiling((double)n * Ni[j] / (double)N);// ni++;
                for (int i = start; i < start + ni && i<n; i++)
                    y[i] =r[rd.Next(start2,start2+Ni[j]-1),1];
                start += ni;
                start2 += Ni[j];
            }
        }

        public static void DoubleSampleNotP(int[,] r, int N, int n1, int n2,out List<int> x2, out int[] y)
        {
            int temp0, temp1;
            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            List<int> Ni = new List<int>();
            List<int> Ind = new List<int>();
            int k = 1;
            Ind.Add(0);
            for (int i = 0; i < N - 1; i++)
            {
                if (r[i, 0] == r[i + 1, 0]) k++;
                else
                {
                    Ni.Add(k); k = 1; Ind.Add(i + 1);
                }
            }


            // 1-й этап, отбор по школам
            List<int> x = new List<int>();
            for (int i = 0; i < n1; i++)
                x.Add(i+1);
            int[] x1 = new int[n1];
            RandomSampleNotP(x, n1, x1);
            

            x2 = new List<int>();
            List<int> x3 = new List<int>();
            for (int i = 0; i < n1; i++)
                for (int j = 0; j < Ni[x1[i]]; j++)
                {
                    x2.Add(Ind[x1[i]] + j);
                    x3.Add(Ind[x1[i]] + j);
                }
            // 2-й этап
            Random rd = new Random();
            if (x3.Count() < n2) y = new int[x3.Count()]; else y = new int[n2];
            int k1;
            for (int i = 0; i < n2&& x3.Count>0; i++)
            {
                k1 = rd.Next(x3.Count());
                y[i] = r[x3[k1], 1];
                x3.RemoveAt(k1);
            }

        }

        public static void DoubleSampleP(int[,] r, int N, int n1, int n2,out List<int> x2 , out int[] y)
        {
            int temp0, temp1;
            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            List<int> Ni = new List<int>();
            List<int> Ind = new List<int>();
            int k = 1;
            Ind.Add(0);
            for (int i = 0; i < N - 1; i++)
            {
                if (r[i, 0] == r[i + 1, 0]) k++;
                else
                {
                    Ni.Add(k); k = 1; Ind.Add(i+1);
                }
            }


            // 1-й этап, отбор по школам
            int[] x1 = new int[n1];
            RandomSampleP(Ni.Count(), n1, x1);

           x2 = new List<int>();
            for (int i = 0; i < n1; i++)
                for (int j = 0; j < Ni[x1[i]]; j++)
                    x2.Add(Ind[x1[i]]+j);
            // 2-й этап
            Random rd = new Random();
            y = new int[n2];
            for (int i = 0; i < n2; i++)
                y[i] = r[x2[rd.Next(n1)],1];
        }

        public static void SerialSampleP(int SN, int Sn, List<int> x, out int[] y, out List<int> index)
        {
            index = new List<int>();
            int[] s=new int[Sn];// номера серий
            RandomSampleP(SN, Sn, s);
            int k = x.Count() / SN; // кол-во элементов в серии
            y = new int[k*Sn];
            for (int i = 0; i < Sn; i++)
            {
                index.Add(x[s[i] * k]);
                if ((s[i] - 1) * k + k < x.Count()) index.Add(x[(s[i] - 1) * k + k - 1]); else index.Add(x[x.Count() - 1]);
                for (int j = 0; j < k && ((s[i] - 1) * k + j < x.Count()); j++)
                    y[i * k + j] = x[(s[i] - 1) * k + j];
            }
        }


        public static void SerialSampleNotP(int SN, int Sn, List<int> x, out int[] y, out List<int> index)
        {
            index = new List<int>();
            int[] s = new int[Sn];// номера серий
            List<int> xs=new List<int>();
            for(int i=1; i<=SN;i++)
                xs.Add(i);
            RandomSampleNotP(xs, Sn, s);
            int k = x.Count() / SN; // кол-во элементов в серии
            y = new int[k * Sn];
            for (int i = 0; i < Sn; i++)
            {
                index.Add(x[(s[i]-1) * k]);
                if ((s[i] - 1) * k + k < x.Count()) index.Add(x[(s[i] - 1) * k + k - 1]); else index.Add(x[x.Count() - 1]);
                for (int j = 0; j < k && ((s[i] - 1) * k + j < x.Count()); j++)
                    y[i * k + j] = x[(s[i] - 1) * k + j];
            }
        }


        public static void TypicalSampleNotP(int N, int[,] r, int n, int[] y)
        {
            int temp0, temp1;
            for (int i = 0; i < N - 1; i++)
                for (int j = i + 1; j < N; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            List<int> Ni = new List<int>();
            int k = 1;
            for (int i = 0; i < N - 1; i++)
            {
                if (r[i, 0] == r[i + 1, 0]) k++;
                else
                {
                    Ni.Add(k); k = 1;
                }
            }


            int ni, start = 0, start2 = 0;
            Random rd = new Random();
            List<int> t = new List<int>();
            for (int j = 0; j < Ni.Count(); j++)
            {
                t.Clear();
                for (int i = start; i < Ni[j]+start; i++)
                    t.Add(r[i,1]);
                start += Ni[j];

                ni = (int)Math.Ceiling((double)n * Ni[j] / (double)N);

                for (int i = start2; i < start2 + ni && i < n; i++)
                {
                    k = rd.Next(t.Count());
                    y[i] = t[k];
                    t.RemoveAt(k);
                }
                start2 += ni;
             }
        }
        public static void MechanicalSample(int N, int n, int[] y)
        {
            int h = N / n;
            for (int i = 0; i < n; i++)
                y[i] = 1+i*h;
        }

        public static void RandomSampleP(int N, int n,int[] y)
        {

            Random rd=new Random();
            for (int i = 0; i < n; i++)
                y[i] = rd.Next(N);
        }
        public static void RandomSampleNotP(List<int> x, int n, int[] y)
        {
            Random rd = new Random();
            int k;
            for (int i = 0; i < n; i++)
            {
                k = rd.Next(x.Count());
                y[i] = x[k];
                x.RemoveAt(k);
            }
        }

        public static double Srednee(int[] x)
        {
            double sr = 0;
            for (int i = 0; i < x.Length; i++)
                sr += x[i];
            return sr/x.Length;
        }
        public static void SredneeTypical(int[,] r,int n,out List<double> sr, out List<int> ind)
        {

            int temp0, temp1;
            for (int i = 0; i < n - 1; i++)
                for (int j = i + 1; j < n; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            ind = new List<int>();
            int k = 1;
            ind.Add(0);
            for (int i = 0; i < n - 1; i++)
            {
                if (r[i, 0] == r[i + 1, 0]) k++;
                else
                {
                     k = 1; ind.Add(i + 1);
                }
            }
            ind.Add(n);
            sr = new List<double>();
            for (int j = 0; j< ind.Count()-1; j++)
            {
                double ss = 0; int p = 0;
                for (int i = ind[j]; i < ind[j + 1]; i++)
                { ss += r[i, 1]; p++; }
                sr.Add(ss / p);
            }

        }
        public static double OtklTypical(List<double> sr, int[,] r, int n, List<int> ind)
        {
            double sko = 0;
            int temp0, temp1;
            for (int i = 0; i < n - 1; i++)
                for (int j = i + 1; j < n; j++)
                    if (r[i, 0] < r[j, 0])
                    {
                        temp0 = r[j, 0]; temp1 = r[j, 1];
                        r[j, 0] = r[i, 0]; r[j, 1] = r[i, 1];
                        r[i, 0] = temp0; r[i, 1] = temp1;
                    }
            double ss = 0; 
            for (int j = 0; j < ind.Count() - 1; j++)
            {
                int p = 0;
                for (int i = ind[j]; i < ind[j + 1]; i++)
                { ss += Math.Pow(r[i, 1] - sr[j], 2); ; p++; }
            }
            return Math.Sqrt(ss/sr.Count());
        }
        public static double Otkl(int[] x)
        {
            double sr = Srednee(x);
            double sko = 0;
            for (int i = 0; i < x.Length; i++)
                sko += Math.Pow(x[i]-sr,2);
            return Math.Sqrt(sko / x.Length);
        }
        public static double OtklSerial(List<double> x)
        {
            double sr = 0;
            for (int i = 0; i < x.Count(); i++)
                sr += x[i];
            sr=sr/x.Count();

            double sko = 0;
            for (int i = 0; i < x.Count(); i++)
                sko += Math.Pow(x[i] - sr, 2);
            return Math.Sqrt(sko / x.Count());
        }
        public static void SredneeSerial(int[,] x, out List<double> sr, List<int> ind)
        {
            sr = new List<double>();
            for (int j = 0; j < ind.Count(); j+=2)
            {
                double ss = 0; int p = 0;
                for (int i = ind[j]; i < ind[j + 1]; i++)
                { ss += x[i,0]; p++; }
                sr.Add(ss / p);
            }
        }

    }
}
