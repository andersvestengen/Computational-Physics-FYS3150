

{
  relerr(i) = log10(arma::max(relerr));
}

string errvals ="errvals.txt";
ofile.open(errvals);

for (int i = 0; i < ex; i++)
{
  cout << relerr(i) << endl;
  ofile << setprecision(15) << nvals(i) <<" "<< relerr(i) << endl;
}
ofile.close();
