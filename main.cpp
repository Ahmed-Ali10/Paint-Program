#if defined(UNICODE) && !defined(_UNICODE)
    #define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
    #define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stack>
#include <vector>
#include <list>

using namespace std;

#define MAXENTRIES 600

#define SAVE 101
#define LOAD 102
#define CLEARSCREEN 103
#define EXIT 104
#define RED 105
#define BLUE 106
#define GREEN 107
#define BLACK 108
#define BCKWHITE 109
#define CURSOR 110

#define LINE_DDA 1
#define LINE_PARA 2
#define LINE_MIDPOINT 3

#define CIRCLE_DIRECT 4
#define CIRCLE_POLAR 5
#define CIRCLE_IPOLAR 6
#define CIRCLE_MIDPOINT 7
#define CIRCLE_MMIDPOINT 8

#define ELIPSE_DIRECT 9
#define ELIPSE_POLAR 10
#define ELIPSE_MIDPOINT 11

#define CIRCLE_FILLING_CIRCLES 12
#define CIRCLE_FILLING_LINES 13

#define FLOOD_FILL 14
#define FLOOD_FILL_NON 15

#define CONVEX 16
#define NON_CONVEX 17

#define CIRCLE_POINT 20
#define CIRCLE_LINE 21

#define REC_POINT 22
#define REC_LINE 23
#define REC_POLYGON 24

#define SQUARE_POINT 25
#define SQUARE_LINE 26

#define HERMITE 27
#define BEZIER 28
#define CARDINAL 29

#define ONE 30
#define TWO 31
#define THREE 32
#define FOUR 33

// Curves and Cardinal Splines
int Round(double x)
{
	return (int)(x + 0.5);
}
struct Vector{
	double v[2];
	Vector(double x = 0, double y = 0)
	 { v[0] = x; v[1] = y; }
	double& operator[](int i){ return v[i];
	}
};

void DrawHermiteCurveModified(HDC hdc,Vector& p1, Vector& T1, Vector& p2, Vector& T2,int xleft, int ytop,int xright,int ybottom,COLORREF c)
{
	double a0 = p1[0], a1 = T1[0],
		a2 = -3 * p1[0] - 2 * T1[0] + 3 * p2[0] - T2[0],
		a3 = 2 * p1[0] + T1[0] - 2 * p2[0] + T2[0];
	double b0 = p1[1], b1 = T1[1],
		b2 = -3 * p1[1] - 2 * T1[1] + 3 * p2[1] - T2[1],
		b3 = 2 * p1[1] + T1[1] - 2 * p2[1] + T2[1];
	for (double t = 0; t <= 1; t += 0.001)
	{
		double t2 = t*t, t3 = t2*t;
		double x = a0 + a1*t + a2*t2 + a3*t3;
		double y = b0 + b1*t + b2*t2 + b3*t3;
		 if(x>=xleft && x<= xright && y>=ytop && y<=ybottom)
        SetPixel(hdc, Round(x), Round(y), c);
	}
}

void DrawHermiteCurve(HDC hdc,Vector& p1, Vector& T1, Vector& p2, Vector& T2,COLORREF c)
{
	double a0 = p1[0], a1 = T1[0],
		a2 = -3 * p1[0] - 2 * T1[0] + 3 * p2[0] - T2[0],
		a3 = 2 * p1[0] + T1[0] - 2 * p2[0] + T2[0];
	double b0 = p1[1], b1 = T1[1],
		b2 = -3 * p1[1] - 2 * T1[1] + 3 * p2[1] - T2[1],
		b3 = 2 * p1[1] + T1[1] - 2 * p2[1] + T2[1];
	for (double t = 0; t <= 1; t += 0.001)
	{
		double t2 = t*t, t3 = t2*t;
		double x = a0 + a1*t + a2*t2 + a3*t3;
		double y = b0 + b1*t + b2*t2 + b3*t3;
        SetPixel(hdc, Round(x), Round(y), c);
	}
}

void DrawBezierCurve(HDC hdc,Vector& P0,Vector& P1,Vector& P2,Vector& P3,int xleft, int ytop,int xright,int ybottom,COLORREF c) {
      Vector T0(3*(P1.v[0]-P0.v[0]),3*(P1.v[1]-P0.v[1]));
      Vector T1(3*(P3.v[0]-P2.v[0]),3*(P3.v[1]-P2.v[1]));
      DrawHermiteCurveModified(hdc,P0,T0,P3,T1,xleft, ytop,xright,ybottom,c);
}

void DrawCardinalSpline(HDC hdc,Vector P[],int n,double c,COLORREF color) {
    double c1=1-c;
    Vector T0(c1*(P[2].v[0]-P[0].v[0]),c1*(P[2].v[1]-P[0].v[1]));
    for(int i=2;i<n-1;i++){
        SetPixel(hdc,P[i].v[0],P[i].v[1],color);
        Vector T1(c1*(P[i+1].v[0]-P[i-1].v[0]),c1*(P[i+1].v[1]-P[i-1].v[1]));
        DrawHermiteCurve(hdc,P[i-1],T0,P[i],T1,color);
        T0=T1;
    }
}

// FLOOD FILL FILLING
struct Vertex{
    int x,y;

    Vertex(int xc=0, int yc=0){
        x = xc;
        y = yc;
    }
};
void RecursiveFloodFill(HDC hdc,int x,int y,COLORREF Bc,COLORREF Fc){
    COLORREF C = GetPixel(hdc,x,y);
    if (C==Bc||C==Fc)return;
    SetPixel(hdc,x,y,Fc);
    RecursiveFloodFill(hdc,x+1,y,Bc,Fc);
    RecursiveFloodFill(hdc,x-1,y,Bc,Fc);
    RecursiveFloodFill(hdc,x,y+1,Bc,Fc);
    RecursiveFloodFill(hdc,x,y-1,Bc,Fc);
}
void NonRecursiveFloodFill(HDC hdc,int x,int y,COLORREF Bc,COLORREF Fc){
    stack<Vertex> s;
    s.push(Vertex(x,y));
    while (!s.empty()){
        Vertex p = s.top();
        s.pop();
        COLORREF C = GetPixel(hdc,p.x,p.y);
        if (C==Bc||C==Fc)
            continue;
        SetPixel(hdc,p.x,p.y,Fc);
        s.push(Vertex(p.x,p.y-1));
        s.push(Vertex(p.x,p.y+1));
        s.push(Vertex(p.x+1,p.y));
        s.push(Vertex(p.x-1,p.y));
    }
}


// Saving and Loading To-From File
bool HDCToFile(const char* FilePath, HDC Context, RECT Area, uint16_t BitsPerPixel = 24)
{
    uint32_t Width = Area.right - Area.left;
    uint32_t Height = Area.bottom - Area.top;
    BITMAPINFO Info;
    BITMAPFILEHEADER Header;
    memset(&Info, 0, sizeof(Info));
    memset(&Header, 0, sizeof(Header));
    Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    Info.bmiHeader.biWidth = Width;
    Info.bmiHeader.biHeight = Height;
    Info.bmiHeader.biPlanes = 1;
    Info.bmiHeader.biBitCount = BitsPerPixel;
    Info.bmiHeader.biCompression = BI_RGB;
    Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
    Header.bfType = 0x4D42;
    Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    char* Pixels = NULL;
    HDC MemDC = CreateCompatibleDC(Context);
    HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
    DeleteObject(SelectObject(MemDC, Section));
    BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
    DeleteDC(MemDC);
    std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
    if (hFile.is_open())
    {
        hFile.write((char*)&Header, sizeof(Header));
        hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
        hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
        hFile.close();
        DeleteObject(Section);
        return true;
    }
    DeleteObject(Section);
    return false;
}
void load(HWND hWnd, HDC &hdc)
{
    char fileName[12] = "picture.bmp";
    if (fileName == "")
        return ;
    HBITMAP hBitmap;
    hBitmap = (HBITMAP)::LoadImage(NULL, fileName, IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    HDC hLocalDC;
    hLocalDC = CreateCompatibleDC(hdc);
    BITMAP qBitmap;
    int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP),reinterpret_cast<LPVOID>(&qBitmap));
    HBITMAP hOldBmp = (HBITMAP)SelectObject(hLocalDC, hBitmap);
    BOOL qRetBlit = BitBlt(hdc, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight,hLocalDC, 0, 0, SRCCOPY);
    SelectObject (hLocalDC, hOldBmp);
    DeleteDC(hLocalDC);
    DeleteObject(hBitmap);
}
void save(HWND &hWnd)
{
    HDC hdc = GetDC(hWnd);
    char fileName[12] = "picture.bmp";
    if (fileName == "")
        return ;
    int windowWidth ;
    int windowHeight;
    RECT rect;
    if(GetWindowRect(hWnd, &rect))
    {
        windowWidth = rect.right - rect.left;
        windowHeight = rect.bottom - rect.top;
    }
    RECT rect1 = {0, 0, windowWidth, windowHeight};
    HDCToFile(fileName,hdc,rect1);
}



int inside_circle(int xc, int yc, int x, int y, int r){

   return (int)pow((x-xc),2.0) + pow((y-yc), 2.0) - (r*r);
}
void PointClippingCircle(HDC hdc, int x, int y, int xc, int yc, int r, COLORREF color){
    if(inside_circle(xc,yc,x,y,r) < 0){
        SetPixel(hdc,x,y,color);
    }
}

// Line Clipping for Circle
void LineClippingCircle(HDC hdc, int x1, int y1, int x2, int y2, int xc, int yc, int r, COLORREF color){
        double dx = x2 - x1;
        double dy = y2 - y1;
        for(double t = 0; t<1; t+=0.001){
            int x = x1 + (t * dx);
            int y = y1 + (t * dy);
            if(inside_circle(xc,yc,x,y,r) < 0)
                SetPixel(hdc, x, y , color);
        }
}

// Point Clipping
void PointClipping(HDC hdc, int x,int y,int xleft,int ytop,int xright,int ybottom, COLORREF color)
{
    if(x>=xleft && x<= xright && y>=ytop && y<=ybottom)
    SetPixel(hdc,x,y,color);
}

// Line Clipping
union OutCode
{
    unsigned All:4;
    struct
    {
        unsigned left:1,top:1,right:1,bottom:1;
    };
};
OutCode GetOutCode(double x,double y,int xleft,int ytop,int xright,int ybottom)
{
    OutCode out;
    out.All=0;
    if(x<xleft)
        out.left=1;
    else if(x>xright)
        out.right=1;
    if(y<ytop)
        out.top=1;
    else if(y>ybottom)
        out.bottom=1;
    return out;
}
void VIntersect(double xs,double ys,double xe,double ye,int x,double *xi,double *yi)
{
    *xi=x;
    *yi=ys+(x-xs)*(ye-ys)/(xe-xs);
}
void HIntersect(double xs,double ys,double xe,double ye,int y,double *xi,double *yi)
{
    *yi=y;
    *xi=xs+(y-ys)*(xe-xs)/(ye-ys);
}
void CohenSuth(HDC hdc,int xs,int ys,int xe,int ye,int xleft,int ytop,int xright,int ybottom)
{
    double x1=xs,y1=ys,x2=xe,y2=ye;
    OutCode out1=GetOutCode(x1,y1,xleft,ytop,xright,ybottom);
    OutCode out2=GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
    while( (out1.All || out2.All) && !(out1.All & out2.All))
    {
        double xi,yi;
        if(out1.All)
        {
            if(out1.left)
                VIntersect(x1,y1,x2,y2,xleft,&xi,&yi);
            else if(out1.top)
                HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out1.right)
                VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else
                HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x1=xi;
            y1=yi;
            out1=GetOutCode(x1,y1,xleft,ytop,xright,ybottom);
        }
        else
        {
            if(out2.left)
                VIntersect(x1,y1,x2,y2,xleft,&xi,&yi);
            else if(out2.top)
                HIntersect(x1,y1,x2,y2,ytop,&xi,&yi);
            else if(out2.right)
                VIntersect(x1,y1,x2,y2,xright,&xi,&yi);
            else
                HIntersect(x1,y1,x2,y2,ybottom,&xi,&yi);
            x2=xi;
            y2=yi;
            out2=GetOutCode(x2,y2,xleft,ytop,xright,ybottom);
        }
    }
    if(!out1.All && !out2.All)
    {

        MoveToEx(hdc,round(x1),round(y1),NULL);
        LineTo(hdc,round(x2),round(y2));
    }
}

// Polygon Clipping
typedef vector<Vertex> VertexList;
typedef bool (*IsInFunc)(Vertex& v,int edge);
typedef Vertex (*IntersectFunc)(Vertex& v1,Vertex& v2,int edge);

VertexList ClipWithEdge(VertexList p,int edge,IsInFunc In,IntersectFunc Intersect)
{
 VertexList OutList;
 Vertex v1=p[p.size()-1];
 bool v1_in=In(v1,edge);
 for(int i=0;i<(int)p.size();i++)
 {
 Vertex v2=p[i];
 bool v2_in=In(v2,edge);
 if(!v1_in && v2_in)
 {
 OutList.push_back(Intersect(v1,v2,edge));
 OutList.push_back(v2);
 }else if(v1_in && v2_in) OutList.push_back(v2);
 else if(v1_in) OutList.push_back(Intersect(v1,v2,edge));
 v1=v2;
 v1_in=v2_in;
 }
return OutList;
}

bool InLeft(Vertex& v,int edge)
{
return v.x>=edge;
}
bool InRight(Vertex& v,int edge)
{
return v.x<=edge;
}
bool InTop(Vertex& v,int edge)
{
return v.y>=edge;
}
bool InBottom(Vertex& v,int edge)
{
return v.y<=edge;
}

Vertex VIntersect(Vertex& v1,Vertex& v2,int xedge)
{
Vertex res;
    res.x=xedge;
    res.y=v1.y+(xedge-v1.x)*(v2.y-v1.y)/(v2.x-v1.x);
    return res;
}
Vertex HIntersect(Vertex& v1,Vertex& v2,int yedge)
{
    Vertex res;
    res.y=yedge;
    res.x=v1.x+(yedge-v1.y)*(v2.x-v1.x)/(v2.y-v1.y);
    return res;
}


void PolygonClip(HDC hdc,POINT *p,int n,int xleft,int ytop,int xright,int ybottom)
{
VertexList vlist;
for(int i=0;i<n;i++)vlist.push_back(Vertex(p[i].x,p[i].y));
    vlist=ClipWithEdge(vlist,xleft,InLeft,VIntersect);
    vlist=ClipWithEdge(vlist,ytop,InTop,HIntersect);
    vlist=ClipWithEdge(vlist,xright,InRight,VIntersect);
    vlist=ClipWithEdge(vlist,ybottom,InBottom,HIntersect);
    Vertex v1=vlist[vlist.size()-1];
for(int i=0;i<(int)vlist.size();i++)
{
    Vertex v2=vlist[i];
    MoveToEx(hdc,round(v1.x),round(v1.y),NULL);
    LineTo(hdc,round(v2.x),round(v2.y));
    v1=v2;
}
}
// CONVEX POLYGON FILLING

struct Entry
{
    int xmin,xmax;
};
void InitEntries(Entry table[])
{
    for(int i=0; i < MAXENTRIES; i++)
        {
            table[i].xmin =  INT_MAX;
            table[i].xmax = - INT_MAX;
        }
}
void ScanEdge(POINT v1,POINT v2,Entry table[])
{
    if(v1.y==v2.y)return;
    if(v1.y>v2.y)
        std::swap(v1,v2);
    double minv=(double)(v2.x-v1.x)/(v2.y-v1.y);
    double x=v1.x;
    int y=v1.y;
    while(y<v2.y){
            if(x<table[y].xmin)table[y].xmin=(int)ceil(x);
            if(x>table[y].xmax)table[y].xmax=(int)floor(x);
            y++;
            x+=minv;
            }
}
void DrawSanLines(HDC hdc,Entry table[],COLORREF color)
{
for(int y=0;y<MAXENTRIES;y++)
if(table[y].xmin<table[y].xmax)
for(int x=table[y].xmin;x<=table[y].xmax;x++)
SetPixel(hdc,x,y,color);
}
void ConvexFill(HDC hdc,POINT p[],int n,COLORREF color)
{
    Entry *table=new Entry[MAXENTRIES];
    InitEntries(table);
    POINT v1=p[n-1];
    for(int i=0;i<n;i++){
            POINT v2=p[i];
            ScanEdge(v1,v2,table);
            v1=p[i];
            }
    DrawSanLines(hdc,table,color);
    delete table;
}

// General Polygon Filling
struct EdgeRec
{
double x;
double minv;
int ymax;
bool operator<(EdgeRec r)
{
return x<r.x;
}
};
typedef list<EdgeRec> EdgeList;
EdgeRec InitEdgeRec(POINT& v1,POINT& v2)
{
    if(v1.y>v2.y)
        std::swap(v1,v2);
    EdgeRec rec;
    rec.x=v1.x;
    rec.ymax=v2.y;
    rec.minv=(double)(v2.x-v1.x)/(v2.y-v1.y);
return rec;
}
void InitEdgeTable(POINT *polygon,int n,EdgeList table[])
{
POINT v1=polygon[n-1];
for(int i=0;i<n;i++)
{
    POINT v2=polygon[i];
    if(v1.y==v2.y){v1=v2;continue;}
    EdgeRec rec=InitEdgeRec(v1, v2);
    table[v1.y].push_back(rec);
    v1=polygon[i];
}
}
void GeneralPolygonFill(HDC hdc,POINT *polygon,int n,COLORREF c)
{

    EdgeList *table=new EdgeList [MAXENTRIES];
    InitEdgeTable(polygon,n,table);
    int y=0;
    while(y<MAXENTRIES && table[y].size()==0)y++;
    if(y==MAXENTRIES)return;
    EdgeList ActiveList=table[y];
    while (ActiveList.size()>0)
        {
        ActiveList.sort();
        for(EdgeList::iterator it=ActiveList.begin();it!=ActiveList.end();it++)
        {
            int x1=(int)ceil(it->x);
            it++;
            int x2=(int)floor(it->x);
            for(int x=x1;x<=x2;x++)SetPixel(hdc,x,y,c);
        }
        y++;
        EdgeList::iterator it=ActiveList.begin();
        while(it!=ActiveList.end())
            if(y==it->ymax) it=ActiveList.erase(it); else it++;
        for(EdgeList::iterator it=ActiveList.begin();it!=ActiveList.end();it++)
            it->x+=it->minv;
        ActiveList.insert(ActiveList.end(),table[y].begin(),table[y].end());
        }
    delete[] table;
}



class Line{
private:
    void swap(int& x1, int& y1, int& x2, int& y2){
        int tmp = x1;
        x1 = x2;
        x2 = tmp;
        tmp = y1;
        y1 = y2;
        y2 = tmp;
    }

    int Round(double x){
        return (int)(x + 0.5);
    }
public:
    Line(){}

    void parametricLine(HDC hdc, int x1, int y1, int x2, int y2,COLORREF c)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        for(double t = 0; t<1; t+=0.001){
            int x = x1 + (t * dx);
            int y = y1 + (t * dy);
            SetPixel(hdc, x, y , c);
        }
    }

    void midPoint(HDC hdc,int x1, int y1, int x2, int y2,COLORREF c)
{


        int dx = x2-x1;
        int dy = y2-y1;

        int x = x1;
        int y = y1;
        int x_dir = 1;
        int y_dir = 1;
        if(dx<0)
        {
            dx = -dx;
            x_dir = -1;
        }
        if(dy<0){
            dy = -dy;
            y_dir = -1;
        }

        if(dy>dx){
            swap(x1,x2,y1,y2);
        }

        if(dx>dy){
            int p =2*dy-dx;
            int two_dy = 2*dy;
            int two_minus= 2*dy - 2*dx;
            for(int i=1; i<=dx; i++){
                x+=x_dir;
                if(p<0)
                    p += two_dy;
                else{
                    p += two_minus;
                    y += y_dir;
                }
                SetPixel(hdc,x,y,c);
            }
        }else{
            int p =2*dx-dy;
            int two_dy = 2*dx;
            int two_minus= 2*dx - 2*dy;
            for(int i=1; i<=dy; i++){
                y+=y_dir;
                if(p<0)
                    p += two_dy;
                else{
                    p += two_minus;
                    x += x_dir;
                }
                SetPixel(hdc,x,y,c);
            }
        }
}
    void DDA(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c){
        int dx = x2 - x1;
        int dy = y2 - y1;
        if (abs(dy) <= abs(dx))
        {
            if (x1 > x2)swap(x1, y1, x2, y2);
            SetPixel(hdc, x1, y1, c);
            int x = x1;
            while (x < x2)
            {
                x++;
                double y = y1 + (double)(x - x1)*dy / dx;
                SetPixel(hdc, x, Round(y), c);
            }
        }
        else {
            if (y1 > y2)swap(x1, y1, x2, y2);
            SetPixel(hdc, x1, y1, c);
            int y = y1;
            while (y < y2)
            {
                y++;
                double x = x1 + (double)(y - y1)*dx / dy;
                SetPixel(hdc, Round(x), y, c);
            }
        }
    }
    ~Line(){}

};

class Circle{

private:
    int Round(double num){
        return (int) num + 0.5;
    }
    void drawPoints(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
        SetPixel(hdc, xc + x, yc + y, c);
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc - x, yc + y, c);
        SetPixel(hdc, xc - x, yc - y, c);
        SetPixel(hdc, xc + y, yc + x, c);
        SetPixel(hdc, xc + y, yc - x, c);
        SetPixel(hdc, xc - y, yc + x, c);
        SetPixel(hdc, xc - y, yc - x, c);
    }
    void drawPointsFourth(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
         SetPixel(hdc, xc + x, yc + y, c);
         SetPixel(hdc, xc + y, yc + x, c);
    }
    void drawPointsThird(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
         SetPixel(hdc, xc - x, yc + y, c);
         SetPixel(hdc, xc -y, yc +x, c);
    }
    void drawPointsSecond(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
         SetPixel(hdc, xc - x, yc - y, c);
         SetPixel(hdc, xc - y, yc - x, c);
    }
    void drawPointsFirst(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc + y, yc - x, c);
    }
    void drawLinesFourth(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
         myLine.parametricLine(hdc, xc, yc ,xc + x, yc + y, c);
         myLine.parametricLine(hdc, xc, yc ,xc + y, yc + x, c);
    }
    void drawLinesThird(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
        myLine.parametricLine(hdc, xc, yc, xc - x ,yc + y, c);
        myLine.parametricLine(hdc, xc, yc, xc - y ,yc + x, c);
    }
    void drawLinesSecond(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
         myLine.parametricLine(hdc,xc,yc, xc - x, yc - y, c);
         myLine.parametricLine(hdc,xc,yc, xc - y, yc - x, c);
    }
    void drawLinesFirst(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
        myLine.parametricLine(hdc,xc,yc, xc + x, yc - y, c);
        myLine.parametricLine(hdc,xc,yc, xc + y, yc - x, c);
    }

public:
    Line myLine;
    Circle(){};

    void Direct(HDC hdc,int xc, int yc, int r, COLORREF c){
        double theta = 0;
        double dtheta = 1.0/r;

        for(theta; theta<360; theta += dtheta){
            double x = xc + r * cos(theta);
            double y = yc + r * sin(theta);

            SetPixel(hdc, Round(x), Round(y), c);
        }
    }

    void Polar(HDC hdc,int xc, int yc, int r, COLORREF c){
        double theta = 0;
        double dtheta = 1.0/r;

        double x = r;
        double y = 0;
        drawPoints(hdc, xc, yc, r, 0, c);

        while(x>y){
            theta += dtheta;
            x = r * cos(theta);
            y = r * sin(theta);
            drawPoints(hdc, xc, yc, Round(x), Round(y), c);
        }
    }

    void Iterative_Polar(HDC hdc,int xc, int yc, int r, COLORREF c){
        double dtheta = 1.0/r;

        double x = r;
        double y = 0;
        double cosinDtheta = cos(dtheta);
        double sinDtheta = sin(dtheta);

        drawPoints(hdc, xc, yc, r, 0, c);

        while(x>y){
            double x1 = x*cosinDtheta - y*sinDtheta;
            double y1 = x*sinDtheta +  y*cosinDtheta;
            x = x1;
            y = y1;
            drawPoints(hdc, xc, yc, Round(x1), Round(y1), c);
        }
    }

    void MidPoint(HDC hdc,int xc, int yc, int r, COLORREF c){
        int x = 0;
        int y = r;
        drawPoints(hdc, xc, yc, r, 0, c);
        while(x<y){
            double d = (x+1)*(x+1) + (y-0.5)*(y-0.5) - r*r;
            if(d>0){
                x++;
                y--;
            }else
            x++;
            drawPoints(hdc, xc, yc, x, y, c);
        }
    }

    void Modified_MidPoint(HDC hdc,int xc, int yc, int r, COLORREF c){
        int x = 0;
        int y = r;
        int d = 1-r;
        drawPoints(hdc, xc, yc, x, y, c);
        while(x<y){
            if(d>0){
                d += 2*(x-y)+5;
                y--;
            }else
                d += 2*x-3;
            x++;
            drawPoints(hdc, xc, yc, x, y, c);
        }
    }

    void FillingWithCircles(HDC hdc, int xc, int yc, int r, int quarter, COLORREF c){
       double theta = 0;
        double dtheta = 1.0/r;

        double x = r;
        double y = 0;

        while(x>y){
            theta += dtheta;
            x = r * cos(theta);
            y = r * sin(theta);
            if(quarter == 1)
                drawPointsFirst(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 2)
                drawPointsSecond(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 3)
                drawPointsThird(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 4)
                drawPointsFourth(hdc, xc, yc, Round(x), Round(y), c);
        }
    }

    void FillingWithLines(HDC hdc, int xc, int yc, int r, int quarter, COLORREF c){
        double x = 0;
        double y = r;
        int rSquared = r*r;
        drawPoints(hdc, xc, yc, 0, r, c);
        while(x<y){
            x++;
            y = sqrt(rSquared - x*x);
            if(quarter == 1)
                drawLinesFirst(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 2)
                drawLinesSecond(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 3)
                drawLinesThird(hdc, xc, yc, Round(x), Round(y), c);
            else if(quarter == 4)
                drawLinesFourth(hdc, xc, yc, Round(x), Round(y), c);
        }
    }

    ~Circle(){};
};

class Elipse{
private:
    int Round(double num){
        return (int) num + 0.5;
    }

    void drawPoints(HDC hdc, int xc, int yc, int x, int y, COLORREF c){
        SetPixel(hdc, xc + x, yc + y, c);
        SetPixel(hdc, xc + x, yc - y, c);
        SetPixel(hdc, xc - x, yc + y, c);
        SetPixel(hdc, xc - x, yc - y, c);
    }

public:
    Elipse(){};
    ~Elipse(){};

    void Direct(HDC hdc,int xc, int yc, int a, int b, COLORREF c){
        double theta = 0;
        double dtheta = a>b?1.0/a:1.0/b;

        for(theta; theta<360; theta += dtheta){
            double x = xc + a * cos(theta);
            double y = yc + b * sin(theta);

            SetPixel(hdc, Round(x), Round(y), c);
        }
    }

    void Polar(HDC hdc,int xc, int yc, int a, int b, COLORREF c){
        double theta = 0;
        double dtheta = a>b?1.0/a:1.0/b;

        double x = a;
        double y = 0;
        drawPoints(hdc, xc, yc, a, 0, c);

        while(x>0 && y<b){
            theta += dtheta;
            x = a * cos(theta);
            y = b * sin(theta);
            drawPoints(hdc, xc, yc, Round(x), Round(y), c);
        }
    }

    void MidPoint(HDC hdc,int xc, int yc, int Rx, int Ry, COLORREF c){
        float P1 = pow(Ry, 2) - pow(Rx, 2)* Ry + (1/4)*pow(Rx,2);
        int x = 0;
        int y = Ry;
        float dx=2*pow(Ry,2)*x;
        float dy=2*pow(Rx,2)*y;
        while (dx < dy)
        {
            if( P1 < 0)
            {
                x++;
                dx=dx+(2 * Ry * Ry);
                P1 += 2 * pow(Ry,2)* x + pow(Ry,2);
            }
            else
            {
                x++;
                y--;
                dx = dx + (2 * Ry * Ry);
                dy = dy - (2 * Rx * Rx);
                P1 += 2 * pow(Ry,2)* x - 2 * pow(Rx, 2)* y + pow(Ry, 2);
            }
            drawPoints(hdc,xc,yc,x,y,c);

        }

        float P2 = pow(Ry, 2) * pow(x + 0.5,2) + pow(Rx, 2) * pow(y-0.5, 2) - pow(Rx,2) * pow(Ry, 2);
        while (y >= 0)
        {
            if (P2 > 0)
            {
                y--;
                dy = dy - (2 * Rx * Rx);
                P2 += pow(Rx,2) - 2*pow(Rx, 2)*y;
            }
            else
            {
                x++;
                y--;
                 dx = dx + (2 * Ry * Ry);
                 dy = dy - (2 * Rx * Rx);
                P2 += 2*pow(Ry, 2)*x - 2*pow(Rx, 2)*y + pow(Rx,2);
            }
            drawPoints(hdc,xc,yc,x,y,c);
        }
    }

};

/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                     HINSTANCE hPrevInstance,
                     LPSTR lpszArgument,
                     int nCmdShow)
{
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_APPSTARTING);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default color as the background of the window */
    wincl.hbrBackground = (HBRUSH) COLOR_BACKGROUND;


    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
           0,                   /* Extended possibilites for variation */
           szClassName,         /* Classname */
           _T("Code::Blocks Template Windows App"),       /* Title Text */
           WS_OVERLAPPEDWINDOW, /* default window */
           CW_USEDEFAULT,       /* Windows decides the position */
           CW_USEDEFAULT,       /* where the window ends up on the screen */
           800,                 /* The programs width */
           500,                 /* and height in pixels */
           HWND_DESKTOP,        /* The window is a child-window to desktop */
           NULL,                /* No menu */
           hThisInstance,       /* Program Instance handler */
           NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}


void showMenus(HWND hwnd){
            HMENU MenuBar = CreateMenu();
            HMENU Line = CreateMenu();
            HMENU Circle = CreateMenu();
            HMENU Elipse = CreateMenu();
            HMENU Filling = CreateMenu();
            HMENU Options = CreateMenu();
            HMENU Color = CreateMenu();
            HMENU CircleFilling = CreateMenu();
            HMENU Clipping = CreateMenu();
            HMENU CircleClipping = CreateMenu();
            HMENU RectangleClipping = CreateMenu();
            HMENU SquareClipping = CreateMenu();
            HMENU Curves = CreateMenu();
            HMENU Quarter = CreateMenu();

            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Options, "Options");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Line, "Line");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Circle, "Circle");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Elipse, "Elipse");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Filling, "Filling");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Clipping, "Clipping");
            AppendMenu(MenuBar, MF_POPUP, (UINT_PTR)Curves, "Curves");

            AppendMenu(Line, MF_STRING, LINE_DDA, "DDA");
            AppendMenu(Line, MF_STRING, LINE_PARA, "Parametric");
            AppendMenu(Line, MF_STRING, LINE_MIDPOINT, "Mid-Point");

            AppendMenu(Circle, MF_STRING, CIRCLE_DIRECT, "Direct");
            AppendMenu(Circle, MF_STRING, CIRCLE_POLAR, "Polar");
            AppendMenu(Circle, MF_STRING, CIRCLE_IPOLAR, "Iterative-Polar");
            AppendMenu(Circle, MF_STRING, CIRCLE_MIDPOINT, "Mid-Point");
            AppendMenu(Circle, MF_STRING, CIRCLE_MMIDPOINT, "Modified Mid-Point");

            AppendMenu(Elipse, MF_STRING, ELIPSE_DIRECT, "Direct");
            AppendMenu(Elipse, MF_STRING, ELIPSE_POLAR, "Polar");
            AppendMenu(Elipse, MF_STRING, ELIPSE_MIDPOINT, "Mid-Point");

            AppendMenu(Options, MF_POPUP, (UINT_PTR)Color, "Choose Color");
            AppendMenu(Options, MF_STRING, BCKWHITE, "Change To White");
            AppendMenu(Options, MF_STRING, CLEARSCREEN, "Clear Screen");
            AppendMenu(Options, MF_STRING, CURSOR, "Change Cursor");
            AppendMenu(Options, MF_STRING, SAVE, "Save");
            AppendMenu(Options, MF_STRING, LOAD, "Load");
            AppendMenu(Options, MF_STRING, EXIT, "Exit");

            AppendMenu(Color, MF_STRING, RED, "Red");
            AppendMenu(Color, MF_STRING, BLUE, "Blue");
            AppendMenu(Color, MF_STRING, GREEN, "Green");
            AppendMenu(Color, MF_STRING, BLACK, "Black");

            AppendMenu(Filling, MF_POPUP, (UINT_PTR)CircleFilling, "Circle");
            AppendMenu(Filling, MF_STRING, FLOOD_FILL, "Flood Fill");
            AppendMenu(Filling, MF_STRING, FLOOD_FILL_NON, "Flood Fill Non Recursive");
            AppendMenu(Filling, MF_STRING, CONVEX, "Convex Filling");
            AppendMenu(Filling, MF_STRING, NON_CONVEX, "Non-Convex Filling");

            AppendMenu(CircleFilling, MF_POPUP, (UINT_PTR)Quarter, "Choose Quarter");
            AppendMenu(CircleFilling, MF_STRING, CIRCLE_FILLING_CIRCLES, "With Circles");
            AppendMenu(CircleFilling, MF_STRING, CIRCLE_FILLING_LINES, "With Lines");

            AppendMenu(Quarter, MF_STRING, ONE, "First Quarter");
            AppendMenu(Quarter, MF_STRING, TWO, "Second Quarter");
            AppendMenu(Quarter, MF_STRING, THREE, "Third Quarter");
            AppendMenu(Quarter, MF_STRING, FOUR, "Fourth Quarter");


            AppendMenu(Clipping, MF_POPUP, (UINT_PTR)CircleClipping, "Circle");
            AppendMenu(Clipping, MF_POPUP, (UINT_PTR)RectangleClipping, "Rectangle");
            AppendMenu(Clipping, MF_POPUP, (UINT_PTR)SquareClipping, "Square");

            AppendMenu(CircleClipping, MF_STRING, CIRCLE_POINT, "Point");
            AppendMenu(CircleClipping, MF_STRING, CIRCLE_LINE, "Line");

            AppendMenu(RectangleClipping, MF_STRING, REC_POINT, "Point");
            AppendMenu(RectangleClipping, MF_STRING, REC_LINE, "Line");
            AppendMenu(RectangleClipping, MF_STRING, REC_POLYGON, "Polygon");

            AppendMenu(SquareClipping, MF_STRING, SQUARE_POINT, "Point");
            AppendMenu(SquareClipping, MF_STRING, SQUARE_LINE, "Line");

            AppendMenu(Curves, MF_STRING, HERMITE, "Square With Hermite");
            AppendMenu(Curves, MF_STRING, BEZIER, "Rectangle With Bezier");
            AppendMenu(Curves, MF_STRING, CARDINAL, "Cardinal Splines");


            SetMenu(hwnd, MenuBar);
}
/*  This function is called by the Windows function DispatchMessage()  */

LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    HDC hdc = GetDC(hwnd);
    static int LineX1,LineX2,LineY1,LineY2;
    static int CircleX1,CircleX2,CircleY1,CircleY2;
    static int ElipseX1,ElipseX2,ElipseX3,ElipseY1,ElipseY2,ElipseY3;
    static int SquareX1,SquareX2,SquareY1,SquareY2;
    static int RecX1,RecX2,RecY1,RecY2,RecX3,RecY3;
    static int LineIndex=0, CircleIndex=0, ElipseIndex=0, SquareIndex=0,RecIndex=0,CurveIndex=0;
    static int PolygonIndex=0, quarter=0;
    static POINT points[5];
    static Vector pointsCurve[4];
    static Vector splinesPoints[6];
    static Line myLine;
    static Circle myCircle;
    static Elipse myElipse;
    static int globalChoice;
    static COLORREF color;
    static boolean CircleDrawn = false;
    static boolean ElipseDrawn = false;
    static HCURSOR hcurs1;
    static HBRUSH hbrush;
    switch (message)                  /* handle the messages */
    {
        case WM_CTLCOLORSTATIC:
            SetBkColor(hdc, RGB(0, 0, 0));
            return (INT_PTR)CreateSolidBrush(RGB(0, 0, 0));
            break;
        case WM_CREATE:
            showMenus(hwnd);
            break;

        case WM_COMMAND:
            switch(LOWORD(wParam)) {

                case RED:
                    color=RGB(255,0,0);
                    break;
                case BLUE:
                    color=RGB(0,0,255);
                    break;
                case GREEN:
                    color=RGB(0,255,0);
                    break;
                case BLACK:
                    color=RGB(0,0,0);
                    break;
                case CLEARSCREEN:
                    char ans;
                    cout << "Are you want to clear window? (y/n)" << endl;
                    cin >> ans;
                    if(ans == 'y')
                        InvalidateRect(hwnd, NULL,TRUE);
                    break;
                case SAVE:
                    save(hwnd);
                    break;
                case LOAD:
                    load(hwnd,hdc);
                    break;
                case EXIT:
                    exit(0);
                    break;
                case ONE:
                    globalChoice=-1;
                    quarter = 1;
                    break;
                case TWO:
                    globalChoice=-1;
                    quarter = 2;
                    break;
                case THREE:
                    globalChoice=-1;
                    quarter = 3;
                    break;
                case FOUR:
                    globalChoice=-1;
                    quarter = 4;
                    break;
                case LINE_DDA:
                    globalChoice = LINE_DDA;
                    LineIndex=0;
                    break;
                case LINE_PARA:
                    globalChoice = LINE_PARA;
                    LineIndex=0;
                    break;
                case LINE_MIDPOINT:
                    globalChoice = LINE_MIDPOINT;
                    LineIndex=0;
                    break;
                case CIRCLE_DIRECT:
                    globalChoice = CIRCLE_DIRECT;
                    CircleIndex=0;
                    break;
                case CIRCLE_POLAR:
                    globalChoice = CIRCLE_POLAR;
                    CircleIndex=0;
                    break;
                case CIRCLE_IPOLAR:
                    globalChoice = CIRCLE_IPOLAR;
                    CircleIndex=0;
                    break;
                case CIRCLE_MIDPOINT:
                    globalChoice = CIRCLE_MIDPOINT;
                    CircleIndex=0;
                    break;
                case CIRCLE_MMIDPOINT:
                    globalChoice = CIRCLE_MMIDPOINT;
                    CircleIndex=0;
                    break;
                case ELIPSE_DIRECT:
                    globalChoice = ELIPSE_DIRECT;
                    ElipseIndex=0;
                    break;
                case ELIPSE_POLAR:
                    globalChoice = ELIPSE_POLAR;
                    ElipseIndex=0;
                    break;
                case ELIPSE_MIDPOINT:
                    globalChoice = ELIPSE_MIDPOINT;
                    ElipseIndex=0;
                    break;
                case CIRCLE_FILLING_CIRCLES:
                    globalChoice = CIRCLE_FILLING_CIRCLES;
                    break;
                case CIRCLE_FILLING_LINES:
                    globalChoice = CIRCLE_FILLING_LINES;
                    break;
                case FLOOD_FILL:
                    globalChoice = FLOOD_FILL;
                    break;
                case FLOOD_FILL_NON:
                    globalChoice = FLOOD_FILL_NON;
                    break;
                case CONVEX:
                    PolygonIndex=0;
                    globalChoice = CONVEX;
                    break;
                case NON_CONVEX:
                    PolygonIndex=0;
                    globalChoice = NON_CONVEX;
                    break;
                case REC_POINT:
                    RecIndex=0;
                    globalChoice = REC_POINT;
                    break;
                case REC_LINE:
                    RecIndex=0;
                    globalChoice = REC_LINE;
                    break;
                case REC_POLYGON:
                    PolygonIndex=0;
                    RecIndex=0;
                    globalChoice = REC_POLYGON;
                    break;
                case SQUARE_POINT:
                    SquareIndex=0;
                    globalChoice = SQUARE_POINT;
                    break;
                case SQUARE_LINE:
                    SquareIndex=0;
                    globalChoice = SQUARE_LINE;
                    break;
                case CIRCLE_POINT:
                    CircleIndex=0;
                    globalChoice = CIRCLE_POINT;
                    break;
                case CIRCLE_LINE:
                    CircleIndex=0;
                    globalChoice = CIRCLE_LINE;
                    break;
                case HERMITE:
                    SquareIndex=0;
                    globalChoice = HERMITE;
                    break;
                case BEZIER:
                    RecIndex=0;
                    globalChoice = BEZIER;
                    break;
                case CARDINAL:
                    CurveIndex=0;
                    globalChoice = CARDINAL;
                    break;
                case BCKWHITE:
                    hbrush = CreateSolidBrush(RGB(255,255,255));
                    SetClassLongPtr(hwnd,GCLP_HBRBACKGROUND,HandleToLong(hbrush));
                    InvalidateRect(hwnd,NULL,TRUE);
                    break;
                case CURSOR:
                    hcurs1 = LoadCursor(NULL, IDC_CROSS);
                    SetClassLongPtr(hwnd,GCLP_HCURSOR,(LONG_PTR) hcurs1);
                    InvalidateRect(hwnd, NULL,TRUE);
                    break;

           }

        case WM_LBUTTONDOWN:
            // LINE ALGORITHMS CASE
            if(globalChoice == 1 || globalChoice == 2 || globalChoice == 3){
                if(LineIndex == 1){
                    LineX1 = LOWORD(lParam);
                    LineY1 = HIWORD(lParam);
                }else if(LineIndex == 2){
                    LineX2 = LOWORD(lParam);
                    LineY2 = HIWORD(lParam);
                    LineIndex=0;
                    if(globalChoice == 1){
                        myLine.DDA(hdc,LineX1,LineY1,LineX2,LineY2,color);
                    } else if(globalChoice == 2){
                        myLine.parametricLine(hdc,LineX1,LineY1,LineX2,LineY2,color);
                    }else if(globalChoice == 3){
                        myLine.midPoint(hdc,LineX1,LineY1,LineX2,LineY2,color);
                    }
                }
            LineIndex++;
            }
            // CIRCLE ALGORITHMS CASE
            else if(globalChoice == 4 || globalChoice == 5 || globalChoice == 6 || globalChoice == 7 || globalChoice == 8){
                if(CircleIndex == 1){
                    CircleX1 = LOWORD(lParam);
                    CircleY1 = HIWORD(lParam);
                }else if(CircleIndex == 2){
                    CircleX2 = LOWORD(lParam);
                    CircleY2 = HIWORD(lParam);
                    CircleDrawn=true;
                    CircleIndex=0;
                    int r = sqrt(pow(CircleX1-CircleX2, 2.0)+pow(CircleY1-CircleY2, 2.0));
                    if(globalChoice == 4){
                        myCircle.Direct(hdc,CircleX1,CircleY1,r,color);
                    } else if(globalChoice == 5){
                        myCircle.Polar(hdc,CircleX1,CircleY1,r,color);
                    }else if(globalChoice == 6){
                        myCircle.Iterative_Polar(hdc,CircleX1,CircleY1,r,color);
                    }else if(globalChoice == 7){
                        myCircle.MidPoint(hdc,CircleX1,CircleY1,r,color);
                    }else if(globalChoice == 8){
                        myCircle.Modified_MidPoint(hdc,CircleX1,CircleY1,r,color);
                    }
                }
                CircleIndex++;
            }
            // ELIPSE ALGORITHMS CASE
            else if(globalChoice == 9 || globalChoice == 10 || globalChoice == 11){
                if(ElipseIndex == 1){
                    ElipseX1 = LOWORD(lParam);
                    ElipseY1 = HIWORD(lParam);
                }else if(ElipseIndex == 2){
                    ElipseX2 = LOWORD(lParam);
                    ElipseY2 = HIWORD(lParam);
                }else if(ElipseIndex == 3){
                    ElipseX3 = LOWORD(lParam);
                    ElipseY3 = HIWORD(lParam);
                    int a = sqrt(pow(ElipseX1-ElipseX2, 2.0)+pow(ElipseY1-ElipseY2, 2.0));
                    int b = sqrt(pow(ElipseX1-ElipseX3, 2.0)+pow(ElipseY1-ElipseY3, 2.0));
                    ElipseIndex=0;
                    ElipseDrawn = true;
                    if(globalChoice == 9){
                        myElipse.Direct(hdc,ElipseX1,ElipseY1,a,b,color);
                    } else if(globalChoice == 10){
                        myElipse.Polar(hdc,ElipseX1,ElipseY1,a,b,color);
                    }else if(globalChoice == 11){
                        myElipse.MidPoint(hdc,ElipseX1,ElipseY1,a,b,color);
                    }
                }
                ElipseIndex++;
            }
            // FILLING CIRCLES WITH LINES AND CIRCLES
            else if((globalChoice == 12 && CircleDrawn) || (globalChoice == 13 && CircleDrawn)){
                int r = sqrt(pow(CircleX1-CircleX2, 2.0)+pow(CircleY1-CircleY2, 2.0));
                if(globalChoice == 12){
                    while(r-->0){
                        myCircle.FillingWithCircles(hdc,CircleX1,CircleY1,r, quarter, color);
                    }
                }
                else if(globalChoice == 13){
                     myCircle.FillingWithLines(hdc, CircleX1,CircleY1,r, quarter, color);
                }
            }
            // FLOOD FILL RECURSIVE CASE
            else if((globalChoice == 14 && CircleDrawn) || (globalChoice == 14 && ElipseDrawn)){
                    if(CircleDrawn)
                        RecursiveFloodFill(hdc,CircleX1, CircleY1, color,color);
                    if(ElipseDrawn)
                        RecursiveFloodFill(hdc,ElipseX1, ElipseY1, color,color);
            }
            // FLOOD FILL NON RECURSIVE CASE
            else if((globalChoice == 15 && CircleDrawn) || (globalChoice == 15 && ElipseDrawn)){
                    if(CircleDrawn)
                        NonRecursiveFloodFill(hdc,CircleX1, CircleY1, color,color);
                    if(ElipseDrawn)
                        NonRecursiveFloodFill(hdc,ElipseX1, ElipseY1, color,color);
            }

            //Convex Filling
            else if(globalChoice == 16){
                 if(PolygonIndex==1)
                {
                   points[0].x = LOWORD(lParam);
                   points[0].y = HIWORD(lParam);
                }
                else if(PolygonIndex==2)
                {
                   points[1].x=LOWORD(lParam);
                   points[1].y=HIWORD(lParam);
                }
                else if(PolygonIndex==3)
                {
                   points[2].x=LOWORD(lParam);
                   points[2].y=HIWORD(lParam);
                }
                else if(PolygonIndex==4)
                {
                   points[3].x=LOWORD(lParam);
                   points[3].y=HIWORD(lParam);
                }
                else if(PolygonIndex==5)
                {
                   points[4].x=LOWORD(lParam);
                   points[4].y=HIWORD(lParam);
                   //change pen color//
                   //SelectObject(hdc, GetStockObject(DC_PEN));
                   //SetDCPenColor(hdc, RGB(60,239,161));
                   Polygon(hdc, points, 5);
                }
                else if(PolygonIndex==6)
                {
                   ConvexFill(hdc, points, 5, color);
                   PolygonIndex=0;
                }
                PolygonIndex++;
            }
            //Non Convex Filling
            else if(globalChoice == 17){
                if(PolygonIndex==1)
                {
                   points[0].x = LOWORD(lParam);
                   points[0].y = HIWORD(lParam);
                }
                else if(PolygonIndex==2)
                {
                   points[1].x=LOWORD(lParam);
                   points[1].y=HIWORD(lParam);
                }
                else if(PolygonIndex==3)
                {
                   points[2].x=LOWORD(lParam);
                   points[2].y=HIWORD(lParam);
                }
                else if(PolygonIndex==4)
                {
                   points[3].x=LOWORD(lParam);
                   points[3].y=HIWORD(lParam);
                }
                else if(PolygonIndex==5)
                {
                   points[4].x=LOWORD(lParam);
                   points[4].y=HIWORD(lParam);
                   //change pen color//
                   //SelectObject(hdc, GetStockObject(DC_PEN));
                   //SetDCPenColor(hdc, RGB(60,239,161));
                   Polygon(hdc, points, 5);
                }
                else if(PolygonIndex==6)
                {
                   GeneralPolygonFill(hdc, points, 5, color);
                   PolygonIndex=0;
                }
                PolygonIndex++;
            }
            // Circle Clipping
            else if(globalChoice == 20 || globalChoice == 21){
                if(CircleIndex == 1){
                    CircleX1 = LOWORD(lParam);
                    CircleY1 = HIWORD(lParam);
                }else if(CircleIndex == 2){
                    CircleX2 = LOWORD(lParam);
                    CircleY2 = HIWORD(lParam);
                    int r = sqrt(pow(CircleX1-CircleX2, 2.0)+pow(CircleY1-CircleY2, 2.0));
                    myCircle.Direct(hdc,CircleX1,CircleY1,r,color);
                }else if(CircleIndex == 3 && globalChoice == 20){
                    int X,Y;
                    CircleIndex = 2;
                    X =  LOWORD(lParam);
                    Y =  HIWORD(lParam);
                    int r = sqrt(pow(CircleX1-CircleX2, 2.0)+pow(CircleY1-CircleY2, 2.0));
                    PointClippingCircle(hdc, X,Y, CircleX1,CircleY1,r,color);
                }else if(CircleIndex == 3 && globalChoice == 21){
                    LineX1 = LOWORD(lParam);
                    LineY1 = HIWORD(lParam);
                }else if(CircleIndex == 4 && globalChoice == 21){
                    LineX2 = LOWORD(lParam);
                    LineY2 = HIWORD(lParam);
                    CircleIndex = 2;
                    int r = sqrt(pow(CircleX1-CircleX2, 2.0)+pow(CircleY1-CircleY2, 2.0));
                    LineClippingCircle(hdc,LineX1,LineY1,LineX2,LineY2,CircleX1,CircleY1,r,color);
                }
                CircleIndex++;
            }
            // Rectangle Clipping
            else if(globalChoice == 22 || globalChoice == 23 || globalChoice == 24){
                 if(RecIndex == 1){
                    RecX1 = LOWORD(lParam);
                    RecY1 = HIWORD(lParam);
                }else if(RecIndex == 2){
                    RecX2 = LOWORD(lParam);
                    RecY2 = HIWORD(lParam);
                }else if(RecIndex == 3){
                    RecX3 = LOWORD(lParam);
                    RecY3 = HIWORD(lParam);
                    RecX2 = RecX1;
                    RecY3 = RecY2;
                    int x4=RecX3, y4=RecY1;
                    PolygonIndex=0;
                    myLine.DDA(hdc, RecX1,RecY1,RecX2,RecY2,color);
                    myLine.DDA(hdc, RecX2,RecY2,RecX3,RecY3,color);
                    myLine.DDA(hdc, RecX3,RecY3,x4,y4,color);
                    myLine.DDA(hdc, x4,y4,RecX1,RecY1,color);
                }
                else if(RecIndex == 4 && globalChoice == 22){
                    RecIndex=3;
                    int X,Y;
                    X = LOWORD(lParam);
                    Y = HIWORD(lParam);
                    PointClipping(hdc,X,Y,RecX3,RecY1,RecX1,RecY2,color);
                }
                else if(RecIndex == 4 && globalChoice == 23){
                    LineX1 = LOWORD(lParam);
                    LineY1 = HIWORD(lParam);
                }
                else if(RecIndex == 5 && globalChoice == 23){
                    LineX2 = LOWORD(lParam);
                    LineY2 = HIWORD(lParam);
                    RecIndex=3;
                    CohenSuth(hdc,LineX1,LineY1,LineX2,LineY2,RecX3,RecY1,RecX1,RecY2);
                }
               else if(PolygonIndex==1 && globalChoice == 24)
                {
                   points[0].x = LOWORD(lParam);
                   points[0].y = HIWORD(lParam);
                }
                else if(PolygonIndex==2&& globalChoice == 24)
                {
                   points[1].x=LOWORD(lParam);
                   points[1].y=HIWORD(lParam);
                }
                else if(PolygonIndex==3&& globalChoice == 24)
                {
                   points[2].x=LOWORD(lParam);
                   points[2].y=HIWORD(lParam);
                }
                else if(PolygonIndex==4&& globalChoice == 24)
                {
                   points[3].x=LOWORD(lParam);
                   points[3].y=HIWORD(lParam);
                }
                else if(PolygonIndex==5&& globalChoice == 24)
                {
                   points[4].x=LOWORD(lParam);
                   points[4].y=HIWORD(lParam);
                   PolygonClip(hdc, points, 5,RecX3,RecY1,RecX1,RecY2);
                }
                RecIndex++;
                PolygonIndex++;
            }
            // Square Clipping
            else if(globalChoice == 25 || globalChoice == 26){

                if(SquareIndex == 1){
                    SquareX1 = LOWORD(lParam);
                    SquareY1 = HIWORD(lParam);
                }else if(SquareIndex == 2){
                    SquareX2 = LOWORD(lParam);
                    SquareY2 = HIWORD(lParam);
                    SquareX2 = SquareX1;
                    int x3=abs(SquareY1-SquareY2)+SquareX1, x4=x3, y3=SquareY1, y4=SquareY2;
                    myLine.DDA(hdc, SquareX1,SquareY1,SquareX2,SquareY2,color);
                    myLine.DDA(hdc, SquareX1,SquareY1,x3,y3,color);
                    myLine.DDA(hdc, SquareX2,SquareY2,x4,y4,color);
                    myLine.DDA(hdc, x3,y3,x4,y4,color);
                }else if(SquareIndex == 3 && globalChoice == 25){
                    SquareIndex = 2;
                    int X,Y;
                    X = LOWORD(lParam);
                    Y = HIWORD(lParam);
                    PointClipping(hdc, X,Y,SquareX1,SquareY1,abs(SquareY1-SquareY2)+SquareX1,SquareY2,color);
                }else if(SquareIndex == 3 && globalChoice == 26){
                    LineX1 = LOWORD(lParam);
                    LineY1 = HIWORD(lParam);
                }else if(SquareIndex == 4 && globalChoice == 26){
                    LineX2 = LOWORD(lParam);
                    LineY2 = HIWORD(lParam);
                    SquareIndex = 2;
                    CohenSuth(hdc,LineX1,LineY1,LineX2,LineY2,SquareX1,SquareY1,abs(SquareY1-SquareY2)+SquareX1,SquareY2);
                }
                SquareIndex++;
            }
            // Filling with Hermite Curves
            else if(globalChoice == 27){
                if(SquareIndex == 1){
                    SquareX1 = LOWORD(lParam);
                    SquareY1 = HIWORD(lParam);
                }else if(SquareIndex == 2){
                    SquareX2 = LOWORD(lParam);
                    SquareY2 = HIWORD(lParam);
                    SquareX2 = SquareX1;
                    SquareIndex=0;
                    int x3=abs(SquareY1-SquareY2)+SquareX1, x4=x3, y3=SquareY1, y4=SquareY2;
                    myLine.DDA(hdc, SquareX1,SquareY1,SquareX2,SquareY2,color);
                    myLine.DDA(hdc, SquareX1,SquareY1,x3,y3,color);
                    myLine.DDA(hdc, SquareX2,SquareY2,x4,y4,color);
                    myLine.DDA(hdc, x3,y3,x4,y4,color);

                    pointsCurve[0].v[0] = SquareX1-20;         pointsCurve[0].v[1] = SquareY1;
                    pointsCurve[1].v[0] = SquareX1-10;         pointsCurve[1].v[1] = (abs(SquareY1-SquareY2)/3)+SquareX1;
                    pointsCurve[2].v[0] = SquareX1-10;         pointsCurve[2].v[1] = (abs(SquareY1-SquareY2)/3)*2+SquareX1;
                    pointsCurve[3].v[0] = SquareX1-20;         pointsCurve[3].v[1] = SquareY2;

                    int start = SquareX1;
                    while(start < x3+20){
                        Vector T3(3 * (pointsCurve[1][0] - pointsCurve[0][0]), 3 * (pointsCurve[1][1] - pointsCurve[0][1]));
                        Vector T4(3 * (pointsCurve[3][0] - pointsCurve[2][0]), 3 * (pointsCurve[3][1] - pointsCurve[2][1]));
                        DrawHermiteCurveModified(hdc, pointsCurve[0], T3, pointsCurve[3], T4,SquareX1,SquareY1,abs(SquareY1-SquareY2)+SquareX1,SquareY2,color);
                        pointsCurve[0].v[0]++;
                        pointsCurve[1].v[0]++;
                        pointsCurve[2].v[0]++;
                        pointsCurve[3].v[0]++;
                        start++;
                    }
                }
                SquareIndex++;
            }
            // Filling with Bezier Curve
            else if(globalChoice == 28){
                if(RecIndex == 1){
                    RecX1 = LOWORD(lParam);
                    RecY1 = HIWORD(lParam);
                }else if(RecIndex == 2){
                    RecX2 = LOWORD(lParam);
                    RecY2 = HIWORD(lParam);
                }else if(RecIndex == 3){
                    RecX3 = LOWORD(lParam);
                    RecY3 = HIWORD(lParam);
                    RecX2 = RecX1;
                    RecY3 = RecY2;
                    int x4=RecX3, y4=RecY1;
                    RecIndex=0;
                    myLine.DDA(hdc, RecX1,RecY1,RecX2,RecY2,color);
                    myLine.DDA(hdc, RecX2,RecY2,RecX3,RecY3,color);
                    myLine.DDA(hdc, RecX3,RecY3,x4,y4,color);
                    myLine.DDA(hdc, x4,y4,RecX1,RecY1,color);

                    pointsCurve[0].v[0] = RecX1;                                pointsCurve[0].v[1] = RecY1-20;
                    pointsCurve[1].v[0] = (abs(RecX1-RecX2)/3) + RecX1;         pointsCurve[1].v[1] = RecY1-10;
                    pointsCurve[2].v[0] = (abs(RecX1-RecX2)/3)*2+RecX1;         pointsCurve[2].v[1] = RecY1-10;
                    pointsCurve[3].v[0] = x4;                                   pointsCurve[3].v[1] = RecY1-20;

                    int start = RecY1;
                    while(start < RecY2+20){
                        DrawBezierCurve(hdc, pointsCurve[0],pointsCurve[1],pointsCurve[2],pointsCurve[3],RecX3,RecY1,RecX1,RecY2,color);
                        pointsCurve[0].v[1]++;
                        pointsCurve[1].v[1]++;
                        pointsCurve[2].v[1]++;
                        pointsCurve[3].v[1]++;
                        start++;
                    }
                }

                RecIndex++;
            }
            // Cardinal Splines
            else if(globalChoice == 29){
                if(CurveIndex == 6){
                    CurveIndex=0;
                    DrawCardinalSpline(hdc,splinesPoints,6,0.2,color);
                }else{
                    splinesPoints[CurveIndex].v[0] = LOWORD(lParam);
                    splinesPoints[CurveIndex].v[1] = HIWORD(lParam);
                }
                CurveIndex++;
            }
            break;
        case WM_DESTROY:
            PostQuitMessage (0);       /* send a WM_QUIT to the message queue */
            break;
        default:                      /* for messages that we don't deal with */
            return DefWindowProc (hwnd, message, wParam, lParam);
    }
    return 0;
}
