#include "mex.h"
#include <bitset>
#include <string>
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


// fp8 加法函数
int fp8Sub(int data1, int data2) {
    data2 ^= (1 << 7);
    unsigned int mant1, mant2, mantNew; // data1尾数，data2尾数，相加之后的尾数
    int sign1, sign2, signNew;          // data1符号位, data2符号位，累加和符号位
    int exp1, exp2, expNew;             // data1阶码，data2阶码，累加和阶码
    int threeBit = 0, fourBit = 0, lastBit = 0;
    int fp8Res;
    unsigned int mantToExp;             // 尾数累加有进位
    int data1_issubnormal = 0, data2_issubnormal = 0; // data1，data2是否是非规格化的数据，0代表不是，1代表是

    // 获取符号位、阶码和尾数
    sign1 = (data1 >> 7) & 0x01;
    sign2 = (data2 >> 7) & 0x01;
    exp1 = ((data1 & 0x7f) >> 3) & 0x0f;
    exp2 = ((data2 & 0x7f) >> 3) & 0x0f;
    mant1 = data1 & 0x07;
    mant2 = data2 & 0x07;

    // 1. 检查操作数中是否有0、Inf、NaN
    if (isNan(data1) || isNan(data2)) {
        return 0x7f;  // NAN 和任何数相加为 NAN    0 1111 111
    }

    if (isZero(data1))
        return data2;                   // 0 + x = x    
    else if (isZero(data2))
        return data1;

    if (isSubnormal(data1)) {
        data1_issubnormal = 1;
        exp1 = 1;  // 非规格化值 指数位为 1-b = 1-7 = -6；
        mant1 = (mant1 | 0x00) << 15;  // 或0 代表没有前导数1
    }
    else
        mant1 = (mant1 | 0x08) << 15;  // 或0x08是加上前导数1，左移15位是为了，低阶向高阶对阶时，尾数会向右移（小数点会左移），这是为了保护有效位不丢失，高低指数若相差14位可忽略不计

    if (isSubnormal(data2)) {
        data2_issubnormal = 1;
        exp2 = 1;
        mant2 = (mant2 | 0x00) << 15;
    }
    else
        mant2 = (mant2 | 0x08) << 15;

    // 3. 进行对阶操作，确定新exp，后续只需要进行尾数的累加
    if (exp1 > exp2) {
        expNew = exp1;
        mant2 >>= (exp1 - exp2);
        exp2 = exp1;
    }
    else {
        expNew = exp2;
        mant1 >>= (exp2 - exp1);
        exp1 = exp2;
    }

    // 4. 尾数累加
    if (sign1 ^ sign2) {  // 一正一负
        if (mant1 > mant2) {
            mantNew = mant1 - mant2;
            signNew = sign1;
        }
        else {
            mantNew = mant2 - mant1;
            signNew = sign2;
        }
    }
    else {  // 同正或同负
        mantNew = mant1 + mant2;
        signNew = sign1;
    }

    // 4. 判断有没有产生进位，并进行规格化
    mantToExp = (mantNew >> 18) & 0x03;  // 尾数本身3位，添加了最高位（或0x08）最开始左移15位，所以右移18位，最后取到最高位

    // 5. 判断有没有产生进位，并进行规格化
    if (mantToExp & 0x02) {     // 有进位
        expNew++;
        lastBit = mantNew & 0x01;
        mantNew >>= 1;
    }
    else if (mantToExp == 0) // 没有前导数1
    {    // 没有进位
        if (mantNew == 0)               // 尾数全为0
        {
            return 0;
        }
        else {
            for (int i = 1; i <= 18; i++)  // 规格化
            {
                expNew--;
                mantNew <<= 1;
                if ((mantNew >> 18) & 0x01)   // 判断前导数是否为1
                    break;
            }  // ok
        }
    }

    // 6. Round to nearest even
    threeBit = (mantNew >> 15) & 0x01;  // 第3位尾数是否为1
    fourBit = (mantNew >> 14) & 0x01;  // 第4位尾数是否为1
    lastBit = (mantNew & ((1 << 14) - 1)) | lastBit; // 后14是否存在1
    mantNew >>= 15;
    if (fourBit && (threeBit || lastBit)) {
        mantNew++;                     // 舍入
    }  // ok

    if (mantNew & 0x10) {             // 舍入之后产生进位
        expNew++;
        mantNew >>= 1;  // ok
    }

    // 7. 溢出判断
    if ((expNew == 15 && mantNew >= 0x0f)||(expNew > 15)) {  // 上溢出
        if (signNew == 1) {
            return 0xfe;  // 返回-448
        }
        else {
            return 0x7e;  // 返回448
        }
    }

    else if (expNew == -1) {  // 非规格化值 
        expNew = 0;
        mantNew >>= 2;
    }

    else if (expNew == 0) {  // 非规格化值 
        expNew = 0;
        mantNew >>= 1;
    }

    else if (expNew == -2) {  // 下溢出
        expNew = 0;
        mantNew >>= 3;
    }
    else if (expNew <= -3) {
        return 0;
    }

    // 构造最终结果
    fp8Res = (((signNew & 0x01) << 7) | ((expNew & 0x0f) << 3) | (mantNew & 0x07));
    return fp8Res;
}

// fp8 加法函数
int fp8Add(int data1, int data2) {
    unsigned int mant1, mant2, mantNew; // data1尾数，data2尾数，相加之后的尾数
    int sign1, sign2, signNew;          // data1符号位, data2符号位，累加和符号位
    int exp1, exp2, expNew;             // data1阶码，data2阶码，累加和阶码
    int threeBit = 0, fourBit = 0, lastBit = 0;
    int fp8Res;
    unsigned int mantToExp;             // 尾数累加有进位
    int data1_issubnormal = 0, data2_issubnormal = 0; // data1，data2是否是非规格化的数据，0代表不是，1代表是

    // 获取符号位、阶码和尾数
    sign1 = (data1 >> 7) & 0x01;
    sign2 = (data2 >> 7) & 0x01;
    exp1 = ((data1 & 0x7f) >> 3) & 0x0f;
    exp2 = ((data2 & 0x7f) >> 3) & 0x0f;
    mant1 = data1 & 0x07;
    mant2 = data2 & 0x07;

    // 1. 检查操作数中是否有0、Inf、NaN
    if (isNan(data1) || isNan(data2)) {
        return 0x7f;  // NAN 和任何数相加为 NAN    0 1111 111
    }

    if (isZero(data1))
        return data2;                   // 0 + x = x    
    else if (isZero(data2))
        return data1;

    if (isSubnormal(data1)) {
        data1_issubnormal = 1;
        exp1 = 1;  // 非规格化值 指数位为 1-b = 1-7 = -6；
        mant1 = (mant1 | 0x00) << 15;  // 或0 代表没有前导数1
    }
    else
        mant1 = (mant1 | 0x08) << 15;  // 或0x08是加上前导数1，左移15位是为了，低阶向高阶对阶时，尾数会向右移（小数点会左移），这是为了保护有效位不丢失，高低指数若相差14位可忽略不计

    if (isSubnormal(data2)) {
        data2_issubnormal = 1;
        exp2 = 1;
        mant2 = (mant2 | 0x00) << 15;
    }
    else
        mant2 = (mant2 | 0x08) << 15;

    // 3. 进行对阶操作，确定新exp，后续只需要进行尾数的累加
    if (exp1 > exp2) {
        expNew = exp1;
        mant2 >>= (exp1 - exp2);
        exp2 = exp1;
    }
    else {
        expNew = exp2;
        mant1 >>= (exp2 - exp1);
        exp1 = exp2;
    }

    // 4. 尾数累加
    if (sign1 ^ sign2) {  // 一正一负
        if (mant1 > mant2) {
            mantNew = mant1 - mant2;
            signNew = sign1;
        }
        else {
            mantNew = mant2 - mant1;
            signNew = sign2;
        }
    }
    else {  // 同正或同负
        mantNew = mant1 + mant2;
        signNew = sign1;
    }

    // 4. 判断有没有产生进位，并进行规格化
    mantToExp = (mantNew >> 18) & 0x03;  // 尾数本身3位，添加了最高位（或0x08）最开始左移15位，所以右移18位，最后取到最高位

    // 5. 判断有没有产生进位，并进行规格化
    if (mantToExp & 0x02) {     // 有进位
        expNew++;
        lastBit = mantNew & 0x01;
        mantNew >>= 1;
    }
    else if (mantToExp == 0) // 没有前导数1
    {    // 没有进位
        if (mantNew == 0)               // 尾数全为0
        {
            return 0;
        }
        else {
            for (int i = 1; i <= 18; i++)  // 规格化
            {
                expNew--;
                mantNew <<= 1;
                if ((mantNew >> 18) & 0x01)   // 判断前导数是否为1
                    break;
            }  // ok
        }
    }

    // 6. Round to nearest even
    threeBit = (mantNew >> 15) & 0x01;  // 第3位尾数是否为1
    fourBit = (mantNew >> 14) & 0x01;  // 第4位尾数是否为1
    lastBit = (mantNew & ((1 << 14) - 1)) | lastBit; // 后14是否存在1
    mantNew >>= 15;
    if (fourBit && (threeBit || lastBit)) {
        mantNew++;                     // 舍入
    }  // ok

    if (mantNew & 0x10) {             // 舍入之后产生进位
        expNew++;
        mantNew >>= 1;  // ok
    }

    // 7. 溢出判断
    if ((expNew == 15 && mantNew >= 0x0f)||(expNew > 15)) {  // 上溢出
        if (signNew == 1) {
            return 0xfe;  // 返回-448
        }
        else {
            return 0x7e;  // 返回448
        }
    }

    else if (expNew == -1) {  // 非规格化值 
        expNew = 0;
        mantNew >>= 2;
    }

    else if (expNew == 0) {  // 非规格化值 
        expNew = 0;
        mantNew >>= 1;
    }

    else if (expNew == -2) {  // 下溢出
        expNew = 0;
        mantNew >>= 3;
    }
    else if (expNew <= -3) {
        return 0;
    }

    // 构造最终结果
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
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("fp8Mul_mex:InvalidInput", "Four input arguments required.");
    }
    
     if (nlhs != 2) {
        mexErrMsgIdAndTxt("fp8Mul_mex:InvalidInput", "Two output arguments required.");
    }
    
    // 从输入参数获取实部和虚部的二进制表示
    int real1 = binaryStringToInt(prhs[0]);
    int imag1 = binaryStringToInt(prhs[1]);
    int real2 = binaryStringToInt(prhs[2]);
    int imag2 = binaryStringToInt(prhs[3]);
    
    // 执行复数乘法
    int temp1 = fp8_Mul(real1,real2);
    int temp2 = fp8_Mul(imag1,imag2);
    int temp3 = fp8_Mul(real1,imag2);
    int temp4 = fp8_Mul(imag1,real2);
    int resultReal = fp8Sub(temp1,temp2);
    int resultImag = fp8Add(temp3,temp4);
 
    // 将结果转换为二进制字符串并返回
    plhs[0] = intToBinaryString(resultReal);
    plhs[1] = intToBinaryString(resultImag);
}


