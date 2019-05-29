template<class Type>
matrix<Type> zero_matrix(int n, int p) {
  matrix<Type> out(n,p);
  out.setZero();
  return out;
}
