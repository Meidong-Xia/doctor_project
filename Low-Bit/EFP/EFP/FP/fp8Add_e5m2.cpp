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

int fp8Add(int data1, int data2) {
    unsigned int mant1, mant2, mantNew; //data1尾数，data2尾数，相加之后的尾数
    int sign1, sign2, signNew;          //data1符号位,data2符号位，累加和符号位
    int exp1, exp2, expNew;             //data1阶码，data2阶码，累加和阶码
    int twoBit = 0, threeBit = 0, lastBit = 0;
    int fp8Res;
    unsigned int mantToExp;				//尾数累加有进位
    int data1_issubnormal = 0, data2_issubnormal = 0; //data1，data2是否是非规格化的数据，0代表不是，1代表是

    // 获取符号位、阶码和尾数
    sign1 = (data1 >> 7) & 0x01;
    sign2 = (data2 >> 7) & 0x01;
    exp1 = ((data1 & 0x7f) >> 2) & 0x1f;
    exp2 = ((data2 & 0x7f) >> 2) & 0x1f;
    mant1 = data1 & 0x03;
    mant2 = data2 & 0x03;

    // 1. 检查操作数中是否有0、Inf、NaN
    if (isNan(data1) || isNan(data2)) {
        return 0x7d;  // NAN 和任何数相加为 NAN    0 111 11 01
    }

    if (isInf(data1) && isInf(data2) && (sign1 ^ sign2))
        return 0x7d;					//Inf + -Inf = NAN
    else if (isInf(data1))
        return data1;					//Inf + x = Inf
    else if (isInf(data2))
        return data2;					//-Inf + x = -Inf
    else if (isZero(data1))
        return data2;					//0 + x = x	
    else if (isZero(data2))
        return data1;

    if (isSubnormal(data1)) {
        data1_issubnormal = 1;
        exp1 = 1;  //非规格化值 指数位为 1-b = 1-15 = -14；
        mant1 = (mant1 | 0x00) << 15;  //或0 代表没有前导数1
    }
    else mant1 = (mant1 | 0x04) << 15;//或0x04是加上前导数1，左移15位是为了，低阶向高阶对阶时，尾数会向右移（小数点会左移），这是为了保护有效位不丢失，高低指数若相差15位可忽略不计

    if (isSubnormal(data2)) {
        data2_issubnormal = 1;
        exp2 = 1;
        mant2 = (mant2 | 0x00) << 15;
    }
    else mant2 = (mant2 | 0x04) << 15;

        
    

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
    if (sign1^sign2) {  // 一正一负
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

    //4.判断有没有产生进位，并进行规格化
    mantToExp = (mantNew >> 17) & 0x03;	//尾数本身2位，添加了最高位（或0x04）最开始左移15位，所以右移17位，最后取到最高位

    // 5. 判断有没有产生进位，并进行规格化
    if (mantToExp & 0x02) {     //有进位
        expNew++;
        lastBit = mantNew & 0x01;
        mantNew >>= 1;
    }
    else if (mantToExp == 0) //没有前导数1
    {    //没有进位
        if (mantNew == 0)				//尾数全为0
        {
            return 0;
        }
        else
        {
            for (int i = 1; i <= 17; i++)//规格化
            {
                expNew--;
                mantNew <<= 1;
                if ((mantNew >> 17) & 0x01)   //判断前导数是否为1
                    break;
            }//ok
        }
    }

    // 6. Round to nearest even
    twoBit = (mantNew >> 15) & 0x01;	//第2位尾数是否为1
    threeBit = (mantNew >> 14) & 0x01;	//第3位尾数是否为1
    lastBit = (mantNew & ((1 << 14) - 1)) | lastBit; //后14是否存在1
    mantNew >>= 15;
    if (threeBit && (twoBit || lastBit)) {
        mantNew++;						//舍入
    }//ok

    if (mantNew & 0x08) {				//舍入之后产生进位
        expNew++;
        mantNew >>= 1;	//ok
    }

    // 7. 溢出判断
    if (expNew >= 31) {  // 上溢出
        if (signNew == 1) {
            return 0xfc;  // 返回负无穷
        }
        else {
            return 0x7c;  // 返回正无穷
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

    else if  (expNew <= -2) {  // 下溢出
       return 0; 
    }

    // 构造最终结果
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
    int result = fp8Add(data1, data2);

    plhs[0] = intToBinaryString(result);
}


