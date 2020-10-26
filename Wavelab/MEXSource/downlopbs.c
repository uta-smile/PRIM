void downlopbs(x,n,lpf,m,y)
double x[],lpf[],y[];
int   n,m;
{
  int n2, i, h, j, k; 
  double s;

  /* lowpass version */
  n2 = n/2;
  k = (m+1)/2 + n - m;
  for(i=0;i<n2;i++){
    s = 0;
    j = k;
    for( h=0; h < m; h++){
      while(j >= n) j -= n;
      while(j < 0) j += n;
      s += lpf[h]*x[j];
      j++;
    }
    k += 2;
    y[i] = s;
  }
}
