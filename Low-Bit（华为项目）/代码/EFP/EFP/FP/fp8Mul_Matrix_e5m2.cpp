#include "mex.h"
#include <bitset>
#include <vector>
#include <string>

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

int fp8Sub(int data1, int data2) {
    data2 ^= (1 << 7);
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

// 将 fp8 字符串转换为整数
int fp8StringToInt(const std::string& fp8String) {
    std::bitset<8> bits(fp8String);
    return static_cast<int>(bits.to_ulong());
}

// 将整数转换为 fp8 字符串
std::string intToFP8String(int integerValue) {
    std::bitset<8> bits(integerValue);
    return bits.to_string();
}

// 从 cell 数组中提取 fp8 矩阵
std::vector<std::vector<int>> extractFP8Matrix(const mxArray* cellArray) {
    mwSize numRows = mxGetM(cellArray);
    mwSize numCols = mxGetN(cellArray);
    std::vector<std::vector<int>> fp8Matrix(numRows, std::vector<int>(numCols, 0));
    
    for (mwIndex i = 0; i < numRows; ++i) {
        for (mwIndex j = 0; j < numCols; ++j) {
            mxArray* cellElement = mxGetCell(cellArray, i + numRows * j);
            std::string fp8String(mxArrayToString(cellElement));
            fp8Matrix[i][j] = fp8StringToInt(fp8String);
        }
    }
    return fp8Matrix;
}

// 从 fp8 矩阵生成 cell 数组
mxArray* createCellArrayFromFP8Matrix(const std::vector<std::vector<int>>& fp8Matrix) {
    mwSize numRows = static_cast<mwSize>(fp8Matrix.size());
    mwSize numCols = static_cast<mwSize>(fp8Matrix[0].size());

    mxArray* cellArray = mxCreateCellMatrix(numRows, numCols);

    for (mwIndex i = 0; i < numRows; ++i) {
        for (mwIndex j = 0; j < numCols; ++j) {
            mxSetCell(cellArray, i + numRows * j, mxCreateString(intToFP8String(fp8Matrix[i][j]).c_str()));
        }
    }

    return cellArray;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 4 || nlhs != 2) {
        mexErrMsgIdAndTxt("fp8Mul_mex:InvalidInput", "Four input and two output arguments required.");
    }

    // 获取输入参数
    mxArray *real1Array = prhs[0];
    mxArray *imag1Array = prhs[1];
    mxArray *real2Array = prhs[2];
    mxArray *imag2Array = prhs[3];

    // 获取矩阵的大小
    mwSize m = mxGetM(prhs[0]);  // Number of rows of A
    mwSize n = mxGetN(prhs[0]);  // Number of columns of A
    mwSize k = mxGetN(prhs[2]);  // Number of columns of B
    
    // 提取 fp8 实部和虚部矩阵  string转int
    std::vector<std::vector<int>> real1 = extractFP8Matrix(prhs[0]);
    std::vector<std::vector<int>> imag1 = extractFP8Matrix(prhs[1]);
    std::vector<std::vector<int>> real2 = extractFP8Matrix(prhs[2]);
    std::vector<std::vector<int>> imag2 = extractFP8Matrix(prhs[3]);

    // 计算乘法并得到结果矩阵
    std::vector<std::vector<int>> resultReal(std::vector<int>m, std::vector<int>k, 0));
    std::vector<std::vector<int>> resultImag(std::vector<int>m, std::vector<int>k, 0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            float sum_real = 0;
            float sum_imag = 0;
            for (int l = 0; l < n; l++) {
                float product_real = half_float::operator*(a_real,b_real)-half_float::operator*(a_imag,b_imag);
                float product_imag = half_float::operator*(a_real,b_imag)+half_float::operator*(a_imag,b_real);
                int temp1 = fp8_Mul(real1[i][j],real2[i][j]);
                int temp2 = fp8_Mul(imag1[i][j],imag2[i][j]);
                int temp3 = fp8_Mul(real1[i][j],imag2[i][j]);
                int temp4 = fp8_Mul(imag1[i][j],real2[i][j]);
                resultReal[i][j] = fp8Sub(temp1,temp2);
                resultImag[i][j] = fp8Add(temp3,temp4);
            }
        }
    }

    // 从结果矩阵生成 cell 数组
    mxArray* resultRealCell = createCellArrayFromFP8Matrix(resultReal);
    mxArray* resultImagCell = createCellArrayFromFP8Matrix(resultImag);

    // 设置输出参数
    plhs[0] = resultRealCell;
    plhs[1] = resultImagCell;
}
