#include "mex.h"
#include <bitset>
using namespace std;

int isNan(int data) {
    int exp = ((data & 0x7f) >> 2) & 0x1f;
    int mant = data & 0x03;
    return ((exp == 0x1f) && (mant != 0)) ? 1 : 0;
}

int isInf(int data) {
    int exp = ((data & 0x7f) >> 2) & 0x1f;
    int mant = data & 0x03;
    return ((exp == 0x1f) && (mant == 0)) ? 1 : 0;
}

int isZero(int data) {
    int exp = ((data & 0x7f) >> 2) & 0x1f;
    int mant = data & 0x03;
    return ((exp == 0) && (mant == 0)) ? 1 : 0;
}

int isSubnormal(int data) {
    int exp = ((data & 0x7f) >> 2) & 0x1f;
    int mant = data & 0x03;
    return ((exp == 0) && (mant != 0)) ? 1 : 0;
}

int fp8_Mul(int data1, int data2)
{
	unsigned int mant1, mant2, mantNew;	//data1尾数，data2尾数，相乘之后的尾数
	unsigned int mantToExp;				//尾数累加有进位
	int sign1, sign2, signNew;			//data1符号位,data2符号位，乘积符号位
	int exp1, exp2, expNew;				//data1阶码，data2阶码，乘积阶码
	int data1_issubnormal = 0, data2_issubnormal = 0; //data1，data2是否是非规格化的数据，0代表不是，1代表是

	int twoBit = 0, threeBit = 0, lastBit = 0;
	int fp8Res;

	// 获取符号位、阶码和尾数
	sign1 = (data1 >> 7) & 0x01;
	sign2 = (data2 >> 7) & 0x01;
	exp1 = ((data1 & 0x7f) >> 2) & 0x1f;
	exp2 = ((data2 & 0x7f) >> 2) & 0x1f;
	mant1 = data1 & 0x03;
	mant2 = data2 & 0x03;

	signNew = sign1 ^ sign2;

	//1. 检查操作数中是否有0、Inf、NaN  没问题
	if (isNan(data1) || isNan(data2)) {
		return 0x7d;					//NAN和任何数相乘为NAN
	}
	if (isInf(data1)) {
		if (isZero(data2)) {				//(+/-)Inf * 0 = NAN
			return 0x7d;
		}
		else if (signNew == 1) {			//Inf * a
			return 0xfc;				//-Inf  1111 1100   
		}
		else if (signNew == 0) {
			return 0x7c;				//+Inf   0111 1100 
		}
	}
	if (isInf(data2)) {
		if (isZero(data1)) {				//(+/-)Inf * 0 = NAN
			return 0x7d;
		}
		else if (signNew == 1) {			//Inf * a
			return 0xfc;				//-Inf
		}
		else if (signNew == 0) {
			return 0x7c;				//+Inf
		}
	}
	if (isZero(data1) || isZero(data2))
		return 0;

	//判断是否为非规格化值	
	if (isSubnormal(data1)) {
		data1_issubnormal = 1;
		exp1 = 1;				 // 非规格化值 指数位为 1-b = 1-15 = -14；
		mant1 = (mant1 | 0x00);  // 或0 代表没有前导数1
	}
	else mant1 = (mant1 | 0x04); // 或0x04是加上前导数1

	if (isSubnormal(data2)) {
		data2_issubnormal = 1;
		exp2 = 1;
		mant2 = (mant2 | 0x00);
	}
	else mant2 = (mant2 | 0x04);

	//2.计算mantNew
	mantNew = mant1 * mant2;

	//3.计算expNew
	expNew = (exp1 - 15) + (exp2 - 15) + 15;




	//判断是否产生进位                          //3位乘3位最多产生6位   0010 0000
	if (mantNew & 0x20) {				//3位乘3位最多产生6位   0010 0000
		expNew++;						//阶码+1
		lastBit = mantNew & 0x01;		//最后一位是否为1
		mantNew >>= 1;					//尾数右移
	}
	else if ((mantNew & 0x10)==0) //没有前导数1
	{    //没有进位
		if (mantNew == 0)				//尾数全为0
		{
			return 0;
		}
		else
		{
			for (int i = 1; i <= 5; i++)//规格化
			{
				expNew--;
				mantNew <<= 1;
				if (mantNew & 0x10)   //判断前导数是否为1
					break;
			}//ok
		}
	}
	//最后5位小数舍去后2位，保留前3位（其中包括1位前导数1）
	//4.Round to nearest even，因为浮点数无法精确表示所有数值，所以要进行舍入处理。
	twoBit = (mantNew >> 2) & 0x01;	//第2位尾数是否为1（从左往右数）
	threeBit = (mantNew >> 1) & 0x01;	//第3位尾数是否为1
	lastBit = (mantNew & 0x01)| lastBit; //3位之后也就是后1位是否存在1
	/*舍入（以保留2位小数为例帮助理解下部分舍入代码，
	两种情况入：a). 0.0011或0.00100001，舍入到最近的数，都得0.01；b). 0.0110两边一样近，向偶数舍入得0.1。
	两种情况舍：a). 0.010001舍入到最近的数，得0.01；b). 0.0010两边一样近，向偶数舍入得0.00）
	*/
    int flag = 0;//对于exp=-2的特殊情况
	if (threeBit && (twoBit || lastBit)) {
		mantNew += 0x04;				//舍入，有效的前3位尾数，+1    0100 0(8)	
	}
	if (mantNew & 0x20) {				//舍入之后产生进位
		expNew++;
		mantNew >>= 1;
        flag = 1;
	}

	if (expNew == -2) {  // 下溢出
		if (mantNew > 16&&flag==0) {
			expNew = 0;
			mantNew = 5;//0100
		}
		else return 0;
	}
	else if (expNew == -1) {  // 非规格化值 
		if (mantNew >= 24) {
			expNew = 0;
			mantNew = 8;//1000
		}
		else {
			expNew = 0;
			mantNew >>= 2;
		}
	}
	else if (expNew == 0) {  // 非规格化值 
		if (mantNew >= 28) {
			expNew = 1;
			mantNew = 0;
		}
        else if (mantNew > 20)
		{
			mantNew = 12;//1100
		}
		else {
			expNew = 0;
			mantNew >>= 1;
		}
	}
	mantNew >>= 2;						//舍去5位的后2位
	//溢出判断
	if (expNew >= 31) {					//上溢出
		if (signNew == 1) {
			return 0xfc;				//返回负无穷
		}
		else {
			return 0x7c;				//返回正无穷
		}
	}
	else if (expNew < 0) {
		return 0;
	}

	fp8Res = (((signNew & 0x01) << 7) | ((expNew & 0x1f) << 2) | (mantNew & 0x03));
	return fp8Res;
}


int binaryStringToInt(const mxArray *binaryString) {
    char *str = mxArrayToString(binaryString);
    std::bitset<8> bits(str); // Assuming 8 bits for this example, adjust as needed
    mxFree(str); // Free the dynamically allocated memory from mxArrayToString
    return static_cast<int>(bits.to_ulong());
}

mxArray *intToBinaryString(int integerValue) {
    std::bitset<8> bits(integerValue); // Assuming 8 bits for this example, adjust as needed
    std::string binaryString = bits.to_string();
    return mxCreateString(binaryString.c_str());
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("fp8Add_mex:InvalidInput", "Two input arguments required.");
    }
    int data1 = binaryStringToInt(prhs[0]);
    int data2 = binaryStringToInt(prhs[1]);

    int result = fp8_Mul(data1, data2);

    plhs[0] = intToBinaryString(result);
}


