#include "mex.h"
#include <bitset>
using namespace std;

// 判断是否为无效数：阶码全1，尾数非全0，则为无效数
int isNan(int data) {
	int exp = ((data & 0x7f) >> 3) & 0x0f;  // 取指数位
	int mant = data & 0x07;                 // 取小数位
	return ((exp == 0x0f) && (mant == 0x07)) ? 1 : 0;
}


// 判断是否为0：阶码全0，尾数全0，则为0，符号位为0，则为正零，符号位为1，则为负零
int isZero(int data) {
	int exp = ((data & 0x7f) >> 3) & 0x0f;  // 取指数位
	int mant = data & 0x07;                 // 取小数位
	return ((exp == 0) && (mant == 0)) ? 1 : 0;
}

// 判断是否为非规格化值：指数位全0，尾数不为0
int isSubnormal(int data) {
	int exp = ((data & 0x7f) >> 3) & 0x0f;  // 取指数位
	int mant = data & 0x07;                 // 取小数位
	return ((exp == 0) && (mant != 0)) ? 1 : 0;
}

int fp8_Mul(int data1, int data2)
{
	unsigned int mant1, mant2, mantNew;	//data1尾数，data2尾数，相乘之后的尾数
	unsigned int mantToExp;				//尾数累加有进位
	int sign1, sign2, signNew;			//data1符号位,data2符号位，乘积符号位
	int exp1, exp2, expNew;				//data1阶码，data2阶码，乘积阶码
	int data1_issubnormal = 0, data2_issubnormal = 0; //data1，data2是否是非规格化的数据，0代表不是，1代表是

	int threeBit = 0, fourBit = 0, lastBit = 0;
	int fp8Res;

	// 获取符号位、阶码和尾数
	sign1 = (data1 >> 7) & 0x01;
	sign2 = (data2 >> 7) & 0x01;
	exp1 = ((data1 & 0x7f) >> 3) & 0x0f;
	exp2 = ((data2 & 0x7f) >> 3) & 0x0f;
	mant1 = data1 & 0x07;
	mant2 = data2 & 0x07;

	signNew = sign1 ^ sign2;

	//1. 检查操作数中是否有0、Inf、NaN  没问题
	if (isNan(data1) || isNan(data2)) {
		return 0x7f;					//NAN和任何数相乘为NAN
	}

	if (isZero(data1) || isZero(data2))
		return 0;

	//判断是否为非规格化值
	if (isSubnormal(data1)) {
		data1_issubnormal = 1;
		exp1 = 1;				 // 非规格化值 指数位为 1-b = 1-15 = -14；
		mant1 = (mant1 | 0x00);  // 或0 代表没有前导数1
	}
	else mant1 = (mant1 | 0x08); // 或0x04是加上前导数1

	if (isSubnormal(data2)) {
		data2_issubnormal = 1;
		exp2 = 1;
		mant2 = (mant2 | 0x00);
	}
	else mant2 = (mant2 | 0x08);

	//2.计算mantNew
	mantNew = mant1 * mant2;

	//3.计算expNew
	expNew = (exp1 - 7) + (exp2 - 7) + 7;




	//判断是否产生进位                          //4位乘4位最多产生8位   1000 0000
	if (mantNew & 0x80) {				//4位乘4位最多产生8位   1000 0000
		expNew++;						//阶码+1
		lastBit = mantNew & 0x01;		//最后一位是否为1
		mantNew >>= 1;					//尾数右移
	}
	else if ((mantNew & 0x40) == 0) //没有前导数1
	{    //没有进位
		if (mantNew == 0)				//尾数全为0
		{
			return 0;
		}
		else
		{
			for (int i = 1; i <= 7; i++)//规格化
			{
				expNew--;
				mantNew <<= 1;
				if (mantNew & 0x40)   //判断前导数是否为1
					break;
			}//ok
		}
	}
	//最后7位小数舍去后3位，保留前4位（其中包括1位前导数1）
	//4.Round to nearest even，因为浮点数无法精确表示所有数值，所以要进行舍入处理。
	threeBit = (mantNew >> 3) & 0x01;	//第3位尾数是否为1（从左往右数）
	fourBit = (mantNew >> 2) & 0x01;	//第4位尾数是否为1
	lastBit = (mantNew & 0x03) | lastBit; //4位之后也就是后2位是否存在1
	//舍入（以保留2位小数为例帮助理解下部分舍入代码，
	//两种情况入：a). 0.0011或0.00100001，舍入到最近的数，都得0.01；b). 0.0110两边一样近，向偶数舍入得0.1。
	//两种情况舍：a). 0.010001舍入到最近的数，得0.01；b). 0.0010两边一样近，向偶数舍入得0.00）

	int flag = 0;//对于exp=-2的特殊情况
	if (fourBit && (threeBit || lastBit)&& expNew > 0) {
		mantNew += 0x08;				//舍入，有效的前4位尾数，+1    1 001 000(8)
	}
	if (mantNew & 0x80) {				//舍入之后产生进位
		expNew++;
		mantNew >>= 1;
		flag = 1;
	}

	if (expNew == -3) {  // 下溢出
		if (mantNew > 64) {  //1 000 000  0.0625
  			expNew = 0;
			mantNew = 8;//1 000
		}
		else
		{
			expNew = 0;
			mantNew >>= 4;
		}
	}
	else if (expNew == -2) {  // 下溢出
		if (mantNew >= 96) {  //1 100 000
			expNew = 0;
			mantNew = 16;//10 000
		}
		else
		{
			expNew = 0;
			mantNew >>= 3;
		}
	}
	else if (expNew == -1) {  // 非规格化值
		if (mantNew >= 0x70) {
			expNew = 0;
			mantNew = 0x20;// 100 000
		}
		else if (mantNew > 0x50) {
			expNew = 0;
			mantNew = 0x18;// 11 000
		}
		else {
			expNew = 0;
			mantNew >>= 2;
		}
	}
	else if (expNew == 0) {  // 非规格化值
		if (mantNew >= 0x78) {
			expNew = 1;
			mantNew = 0;
		}
		else if (mantNew > 0x68)
		{
			mantNew = 0x38;//1100
		}
		else if (mantNew >= 0x58)
		{
			mantNew = 0x30;//1100
		}
		else if (mantNew > 0x48)
		{
			mantNew = 0x28;//1100
		}
		else {
			expNew = 0;
			mantNew >>= 1;
		}
	}
	mantNew >>= 3;						//舍去7位的后3位
	//溢出判断
	if ((expNew == 15 && mantNew >= 0x0f) || (expNew > 15)) {  // 上溢出
		if (signNew == 1) {
			return 0xfe;  // 返回-448
		}
		else {
			return 0x7e;  // 返回448
		}
	}
	else if (expNew < 0) {
		return 0;
	}

	fp8Res = (((signNew & 0x01) << 7) | ((expNew & 0x0f) << 3) | (mantNew & 0x07));
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


