inline double Cn(size_t n,size_t k)
{
	if (!n) throw;
	if(!k||k==n) return 1.0;
	double tmp1=n,tmp2=1;
	for(size_t i=1;i<=k-1;++i)
		tmp1*=(n-i);
	for(size_t i=2;i<=k;++i)
		tmp2*=i;
	return tmp1/tmp2;
}

void Diff(size_t num,std::vector<double> &_rst)  //注意函数没有除以步长a^num
{
	_rst.resize(num+1);
	if(!num) {_rst[0]=1;return;}
	for(size_t i=0;i<num+1;++i)
	{
		if(!((num-i)%2))	_rst[i]=Cn(num,num-i);
		else _rst[i]=-Cn(num,i);
	}
}