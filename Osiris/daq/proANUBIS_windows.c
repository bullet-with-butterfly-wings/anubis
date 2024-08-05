/*******************************************************************************/
/*                                                                             */
/*                          V767 test program                                  */
/*                                                                             */
/*                               C.A.E.N                                       */
/*                                                                             */
/*                 by Carlo Tintori     Viareggio, 07/97                       */
/*                                                                             */    
/*******************************************************************************/


#include <stdio.h>
#include <ctype.h>
#include <math.h> 
#include <conio.h>
#include <Windows.h>
#include <stdbool.h>
#include <sys/timeb.h>
#include <CAENVMElib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <vector>
#include <iomanip>
#include <time.h>
#include <thread>
#include <ctime>

#define ulong  unsigned long
#define ushort unsigned short
#define uchar  unsigned char

#define MAX_BLT_SIZE  (4*1024)

ushort out_buffer = 0x0000;
ushort geo_addr = 0x0004;
ushort bit_set = 0x0006;
ushort bit_clear = 0x0008;
ushort fw_rev = 0x104E;
ushort int_level = 0x000A;
ushort vme_reset = 0x0018;
ushort cont_reg_1 = 0x0010;
ushort stat_reg_1 = 0x000E;
ushort stat_reg_2 = 0x0048;
ushort evt_counter = 0x004C;
ushort op_handshake = 0x0050;
ushort op_reg = 0x0052;
ushort tdc_clear = 0x0054;
ushort soft_trig = 0x005A;
ushort man_id_0 = 0x1026;
ushort man_id_1 = 0x102B;
ushort man_id_2 = 0x102E;
ushort board_id_0 = 0x1032;
ushort board_id_1 = 0x1037;
ushort board_id_2 = 0x103B;
ushort board_id_3 = 0x103E;
ushort mcst_addr = 0x0016;
ushort mcst_ctrl = 0x0020;

//number of TDCs
int nTDCs = 6;
int tdc = 0;

//TDC V767 Base addresses
ulong base_addr[6] = { 0xEE000000, 0xEE010000, 0xEE020000, 0xEE030000, 0xEE040000, 0xEE050000 };

// handle per V767 
ulong handle;
int32_t BHandle;

int VMEerror = 0;
char ErrorString[100];

//Hacky way to stop the program after "deamonizing" it. Reads a cfg file which should only contain a zero or one, stops if it doesn't see a one.
bool checkRunFlag()
{
    std::ifstream flagFile("cfg/runFlag.txt");
    int state = 0;
    flagFile >> state;
    if (state == 1) { return true; }
    return false;
}

/*******************************************************************************/
/*                               READ_BLK                                      */
/*******************************************************************************/
int read_blk(uint32_t* buffer)
{
    int byte_cnt = 0;
    CVErrorCodes ret_blck = CAENVME_BLTReadCycle(handle, base_addr[tdc], (uchar*)buffer, 2048, cvA32_U_MBLT, cvD64, &byte_cnt);
    if (ret_blck != cvSuccess) {
        sprintf(ErrorString, "Block transfer failed for TDC %i.", tdc);
        VMEerror = 1;
    }
    return byte_cnt / 4;
}

/*******************************************************************************/
/*                               READ_REG                                      */
/*******************************************************************************/
ulong read_reg(ushort reg_addr)
{
    int dSize = 16;
    ulong data = 0;
    CVErrorCodes ret;
    if (dSize == 16) {
        ret = CAENVME_ReadCycle(handle, base_addr[tdc] + reg_addr, &data, cvA32_U_DATA, cvD16);
    }
    else {
        ret = CAENVME_ReadCycle(handle, base_addr[tdc] + reg_addr, &data, cvA32_U_DATA, cvD32);
    }
    if (ret != cvSuccess) {
        sprintf(ErrorString, "Cannot read at address %08X\n", (uint32_t)(base_addr[tdc] + reg_addr));
        VMEerror = 1;
    }
    return(data);
}

/*******************************************************************************/
/*                                WRITE_REG                                    */
/*******************************************************************************/
void write_reg(ushort reg_addr, ushort data)
{
    CVErrorCodes ret;
    ret = CAENVME_WriteCycle(handle, base_addr[tdc] + reg_addr, &data, cvA32_U_DATA, cvD16);
    if (ret != cvSuccess) {
        sprintf(ErrorString, "Cannot write at address %08X\n", (uint32_t)(base_addr[tdc] + reg_addr));
        VMEerror = 1;
    }
}

/*******************************************************************************/
/*                                WRITE_OPCODE                                    */
/*******************************************************************************/
int write_opcode(ushort code)
{
    ushort rdata;
    int attempts = 0;
    do
    {
        rdata = read_reg(op_handshake);
        attempts++;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    } while ((rdata != 0x02) && (attempts < 50));
    if (attempts > 50)
    {
        printf("Handshake timeout!\n");
        return -1;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    write_reg(op_reg, code);
    return 0;
}

/*******************************************************************************/
/*                                READ_OPCODE                                    */
/*******************************************************************************/
int read_opcode()
{
    ushort rdata;
    int attempts = 0;
    do
    {
        rdata = read_reg(op_handshake);
        attempts++;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    } while ((rdata != 0x01) && (attempts < 50));
    if (attempts > 50)
    {
        printf("Handshake timeout!\n");
        return -1;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    rdata = read_reg(op_reg);
    return rdata;
}
/*******************************************************************************/
/*                                Interpret_TDCError                                */
/*******************************************************************************/
bool interpret_tdcError(ushort chip)
{
    write_opcode(0x8000 + (chip >> 2));
    ushort errStat = read_opcode();
    return errStat;
}

//Get the current time
std::pair<uint32_t, uint32_t> getTime()
{
    auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch());
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(nanoseconds);
    nanoseconds -= seconds;
    std::pair<uint32_t, uint32_t> outTime;
    outTime.first = seconds.count();
    outTime.second = nanoseconds.count();
    return outTime;
}

/*******************************************************************************/
/*                                READ_TDCError                                */
/*******************************************************************************/
bool read_tdcError()
{
    ushort status = read_reg(stat_reg_2);
    ushort tdcErrorMask = 0x8;
    ushort chip3mask = 0x8000;
    ushort chip2mask = 0x4000;
    ushort chip1mask = 0x2000;
    ushort chip0mask = 0x1000;
    bool tdcError = false;
    if (status & tdcErrorMask) { tdcError = true; }
    if (status & chip0mask)
    {
        tdcError = true;
        printf("TDC Chip A is in error.\n");
    }
    if (status & chip1mask)
    {
        tdcError = true;
        printf("TDC Chip B is in error.\n");
    }
    if (status & chip2mask)
    {
        tdcError = true;
        printf("TDC Chip C is in error.\n");
    }
    if (status & chip3mask)
    {
        tdcError = true;
        printf("TDC Chip D is in error.\n");
    }
    if (!tdcError) { printf("No TDC Chips in Error.\n"); }
    return tdcError;
}

/******************************************************************************/
/*                                   MAIN                                     */
/******************************************************************************/
int main(int argc, char *argv[])
{
    time_t seconds = 7200;
    if (argc > 1) { seconds = atol(argv[1]); }
 
    uint32_t usb_link = 0;
    if (CAENVME_Init(cvV1718, 0, 0, &BHandle) == cvSuccess)
    {
        printf("Connected to the device.\n");
    }
    else {
        printf("\tError in opening V1718 bridge \n");
        return 1;
    }

    if (VMEerror) {
        printf(ErrorString);
        CAENVME_End(handle);
        return 1;
    }

    // ------------------------------------------------------------------------------------
    // Acquisition loop
    // ------------------------------------------------------------------------------------
    printf("Starting Acquisition from V767 modules\n");
    short win_offs = -10;
    short win_width = 50;

    for (; tdc < nTDCs; tdc++) {
        printf("Setting up TDC %i.\n", tdc);
        bool tdcError;
        do {
            if (tdc < 5) //Set up TDCs in triggered mode
            {
                write_reg(vme_reset, 0x0011);
                std::this_thread::sleep_for(std::chrono::milliseconds(2000));
                write_opcode(0x1500); //Load default config
                write_opcode(0xB800); //Set mask for evnt overflow
                write_opcode(0xB);
                write_opcode(0x2400);//Disable all channels
                write_opcode(0x1000);//sets StopTriggerMatching mode
                write_opcode(0x3000);//sets window width
                write_opcode(win_width);
                write_opcode(0x3200);//sets window offset
                write_opcode(win_offs);
                write_opcode(0x7000);//sets Data_Ready_Mode=Event Ready
                write_opcode(0x2300);//Enable all channels
            }
            else //Set up TDC 6 in continuous readout.
            {
                write_reg(vme_reset, 0x0011);
                std::this_thread::sleep_for(std::chrono::milliseconds(2000));
                write_opcode(0x2100);
                //write_opcode(0x2102);
                write_opcode(0x1300);//Set continuous mode
                write_opcode(0x7200);//sets Data_Ready_Mode=Buffer not empty
                //write_opcode(0x1600); //Set this as the default mode
                //write_opcode(0x1800); //Enable auto-load
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            tdcError = read_tdcError();
        } while (tdcError);
    }

    std::pair<uint32_t, uint32_t> thisTime;
    uint32_t thisWord;
    uint32_t fillWord;
    uint32_t buffer[MAX_BLT_SIZE] = { 0 };
    std::vector<uint32_t> dataBuffer;

    while (true) {
        //Save the chrono start time
        thisTime = getTime();

        time_t endwait;
        time_t start = time(NULL);
        time_t shortwaittime = 1800;
        time_t shortwait = start + shortwaittime;
        endwait = start + seconds;
        //enum states {Waiting, Ready, Invalid};
        //states state=Waiting;

        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&now), "%y%m%d_%H%M");
        std::ostringstream date;
        date << std::put_time(std::localtime(&now), "%y_%m");
        std::ostringstream day;
        day << std::put_time(std::localtime(&now), "%y_%m_%d");
        std::string ofName = "data/" + date.str() + "/" + day.str() + "/proAnubis_" + oss.str() + ".raw";
        CreateDirectory(("data/" + date.str()).c_str(), NULL);
        CreateDirectory(("data/" + date.str() + "/" + day.str()).c_str(), NULL);

        std::ofstream thisFile;
        std::ofstream monFile;
        thisFile.open(ofName.c_str(), std::ios::binary | std::ios::out);
        //Write the start time into the TFile
        thisFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t)); 
        thisFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
        monFile.open("tempMonFile.raw", std::ios::binary | std::ios::out);
        monFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
        monFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
        ulong EOBMask = 0x200000;
        ulong headerMask = 0x400000;

        tdc = 0;
        for (; tdc < nTDCs; tdc++)
        {
            write_reg(tdc_clear, 0x0011);
        }
        tdc = 0;

        do {
            int wordCount = read_blk(buffer);
            thisTime = getTime();
            bool lastWasBad = false;
            for (int b = 0; b < wordCount; b++)
            {
                thisWord = buffer[b];
                if (thisWord & headerMask)
                {
                    if (thisWord & EOBMask)
                    {
                        if (lastWasBad) { break; }
                        else
                        {
                            lastWasBad = true;
                            fillWord = thisWord;
                        }
                    }
                    else 
                    {
                        if (lastWasBad)
                        {
                            lastWasBad = false;
                            dataBuffer.push_back(fillWord);
                        }
                        dataBuffer.push_back(thisWord); 
                    }
                }
                else
                {
                    if (lastWasBad)
                    {
                        lastWasBad = false;
                        dataBuffer.push_back(fillWord);
                    }
                    dataBuffer.push_back(thisWord);
                }
            }
            if (dataBuffer.size() > 0)
            {
                thisFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t)); 
                thisFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
                uint32_t headerWord = (tdc << 24) + dataBuffer.size();
                thisFile.write((char*)&headerWord, sizeof(uint32_t));
                thisFile.write((char*)&dataBuffer[0], dataBuffer.size() * sizeof(uint32_t));
                monFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
                monFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
                monFile.write((char*)&headerWord, sizeof(uint32_t));
                monFile.write((char*)&dataBuffer[0], dataBuffer.size() * sizeof(uint32_t));
                dataBuffer.clear();
            }
            tdc++;
            if (tdc >= nTDCs)
            {
                tdc = 0;
                if (!checkRunFlag()) {
                    thisTime = getTime();
                    thisFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t)); //(char*) &a, sizeof(a));
                    thisFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
                    monFile.close();
                    std::remove("data/monitor/proANUBISmonitor.raw");
                    std::rename("tempMonFile.raw", "data/monitor/proANUBISmonitor.raw");
                    thisFile.close();
                    CAENVME_End(handle);
                    return 0;
                }
                if (time(NULL) > shortwait)
                {
                    thisTime = getTime();
                    monFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
                    monFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
                    monFile.close();
                    std::remove("data/monitor/proANUBISmonitor.raw");
                    std::rename("tempMonFile.raw", "data/monitor/proANUBISmonitor.raw");
                    monFile.open("tempMonFile.raw", std::ios::binary | std::ios::out);
                    monFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
                    monFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
                    shortwait = shortwait + shortwaittime;
                }
            }
        } while (time(NULL) < endwait);
        thisTime = getTime();
        thisFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
        thisFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
        thisFile.close();
        monFile.write(reinterpret_cast<char*>(&thisTime.first), sizeof(uint32_t));
        monFile.write(reinterpret_cast<char*>(&thisTime.second), sizeof(uint32_t));
        monFile.close();
        std::remove("data/monitor/proANUBISmonitor.raw");
        std::rename("tempMonFile.raw", "data/monitor/proANUBISmonitor.raw");
    }
    CAENVME_End(handle);
    return 0;
}
