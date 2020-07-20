#include <Xlib.h>
#include <Xatom.h>
#include <stdio.h>
#include <string.h>
#include <Xutil.h>


#define TRUE 1
#define FALSE 0
#define MAXPATH 500
#define PSOnly 0
#define XOnly 1
#define XAndPS 2

extern char **environ;

static FILE *PSFile;

static GC GraphContext, StopContext, StartContext, ContContext,
          StepContext, BackContext; 
static int WinULX, WinULY, LastWinX, LastWinY, CurrentDashLength;
static int PSULX, PSULY, PathLength, CurrentLineWidth, CurrentLineStyle;
static unsigned int WinWidth, WinHeight, BorderWidth, Dep;
static unsigned int PSWidth, PSHeight;
static Window root, StopButton, StartButton, StepButton, ContButton;
static Window GraphWin;
static Display *display;
static XEvent report;
static int screen;
static double UserLLX, UserLLY, UserWidth, UserHeight;
static Colormap Kleurtjes;
static Pixmap Backup;
static Font CurrentFont;
static XFontStruct *CurrentFontStruct;
static XStandardColormap *mapinfo;
static char *CurrentFontName, *PSFileName;
static int CurrentFontSize;
static double Scale;
unsigned long BLUE, GREEN, RED, YELLOW, BLACK, WHITE;
static short int InitOK, FramOK, StartOK, StopOK, Output;
static char *BackColName;



void 
error (Message)
char *Message;
{
   printf ("%s\n", Message);
   exit(1);
}


void 
Framing (Xll, Yll, Xur, Yur)
double Xll, Yll, Xur, Yur;
{
   FramOK = TRUE;
   if ((Xll!=Xur) && (Yll!=Yur)) 
      {
         UserLLX = Xll;
         UserLLY = Yll;
         UserWidth  = Xur-Xll;
         UserHeight = Yur-Yll;
      }
   else
      { printf("ERROR in function FRAMING: invalid arguments .\n ");
        exit(1);
      }
}


void 
MapStartButton()

{
   int Xul, Yul;

   if (StartOK == TRUE) 
      {
         XGetGeometry (display, GraphWin, &root, &WinULX, &WinULY, &WinWidth,
                                         &WinHeight, &BorderWidth, &Dep);
         Xul =  WinULX+0.03*WinWidth;
         Yul =  WinULY+0.06*WinHeight;
         XMoveWindow (display, StartButton, Xul, Yul);
         XMapRaised(display, StartButton);
         XDrawString (display, StartButton, StartContext, 5,12,"start",5); 
      }
}

void 
MapStopButton()

{
   int Xul, Yul;

   if (StopOK == TRUE) 
     {
       XGetGeometry (display, GraphWin, &root, &WinULX, &WinULY, &WinWidth,
                                         &WinHeight, &BorderWidth, &Dep);
       Xul =  WinULX+0.03*WinWidth;
       Yul =  WinULY+0.01*WinHeight;
       XMoveWindow (display, StopButton, Xul, Yul);
       XMapRaised(display, StopButton);
       XDrawString (display, StopButton, StopContext, 5,12,"stop",4); 
/*       Xul =  WinULX+0.03*WinWidth;
       Yul =  WinULY+0.06*WinHeight;
       XMoveWindow (display, StartButton, Xul, Yul);
       XMapRaised(display, StartButton);
       XDrawString (display, StartButton, StartContext, 5,12,"start",5); */
    }
}



void
AssociateColor (ColorName, ColorPixel)
char *ColorName;
unsigned long *ColorPixel;
{
   XColor Kleur, KleurExact;

   XAllocNamedColor (display,XDefaultColormap(display,screen),
                     ColorName, &Kleur, &KleurExact); 
   *ColorPixel = Kleur.pixel;
}

void
InstallFastColors ()
{
   XColor Kleur, KleurExact;

   if (Output == XOnly || Output == XAndPS)
     { 
       AssociateColor ("blue", &BLUE);
       AssociateColor ("green", &GREEN);
       AssociateColor ("red", &RED);
       AssociateColor ("yellow", &YELLOW);
       AssociateColor ("white", &WHITE);
       AssociateColor ("black", &BLACK);
     }
}   
   

void
SetQuickColor(ClPx)
XColor *ClPx;
{
  XAllocColor (display, XDefaultColormap(display,screen), ClPx); 
}


unsigned long
ColorCode(R, G, B)
int R, G, B;
{
  XColor *ColorData;
  unsigned long r;

  ColorData = (XColor *) malloc (sizeof(ColorData));

  ColorData->red   = (unsigned short) (R);
  ColorData->green = (unsigned short) (G);
  ColorData->blue  = (unsigned short) (B);
  ColorData->pixel = 0;
  SetQuickColor(ColorData);
  r = (ColorData->pixel);
  free (ColorData);
  return (r);
}

void
SetNumColor(ClPx)
unsigned long ClPx;
{
  XSetForeground(display, GraphContext, ClPx);
  XSetForeground(display, BackContext, ClPx); 
}

void
SetFastColor (Color)
int Color;

{
   if (Output == XOnly || Output == XAndPS)
     { 
       switch(Color)
       {
         case 1 : XSetForeground (display, GraphContext, BLUE); 
                  XSetForeground (display, BackContext, BLUE); 
                  break;
         case 2 : XSetForeground (display, GraphContext, GREEN); 
                  XSetForeground (display, BackContext, GREEN); 
                  break;
         case 3 : XSetForeground (display, GraphContext, RED); 
                  XSetForeground (display, BackContext, RED); 
                  break;
         case 4 : XSetForeground (display, GraphContext, YELLOW); 
                  XSetForeground (display, BackContext, YELLOW); 
                  break;
         case 5 : XSetForeground (display, GraphContext, WHITE); 
                  XSetForeground (display, BackContext, WHITE); 
                  break;
         case 6 : XSetForeground (display, GraphContext, BLACK); 
                  XSetForeground (display, BackContext, BLACK); 
       }
    }
}


void 
SetNamedColor (ColorName)
char *ColorName;

{
   XColor Kleur, KleurExact;
   unsigned long ColorPix;

   if (Output == XOnly || Output == XAndPS)
     { 
       XAllocNamedColor (display,XDefaultColormap(display,screen),ColorName, &Kleur, &KleurExact); 
       ColorPix = Kleur.pixel;
       XSetForeground (display, GraphContext, ColorPix); 
       XSetForeground (display, BackContext, ColorPix); 
    }
}


void
SetNamedBackground (ColorName)
char *ColorName;

{
   XColor Kleur, KleurExact;
   unsigned long ColorPix;

   if (Output == XOnly || Output == XAndPS)
     { 
       XAllocNamedColor (display,XDefaultColormap(display,screen),ColorName, &Kleur, &KleurExact); 
       ColorPix = Kleur.pixel;
       XSetForeground (display, GraphContext, ColorPix); 
       XFillRectangle (display, GraphWin, GraphContext, 0, 0, WinWidth, WinHeight);
       XSetForeground (display, BackContext, ColorPix); 
       XFillRectangle (display, Backup, BackContext, 0, 0, WinWidth, WinHeight);
       SetNamedColor("black");
    }
}


void
SetFillStyle(style)
char *style;
{
   int FillStyle;

   if (Output == XOnly || Output == XAndPS)
     { 
       if (style[0] == 's' || style[0] == 'S') 
         FillStyle = FillSolid;
       else if (style[0] == 't' || style[0] == 'T') 
         FillStyle = FillTiled;
       else if (style[0] == 'p' || style[0] == 'P') 
         FillStyle = FillStippled;
       else if (style[0] == 'o' || style[0] == 'O') 
         FillStyle = FillOpaqueStippled;
       XSetFillStyle (display, GraphContext, FillStyle);
       XSetFillStyle (display, BackContext, FillStyle);
    }
}


char
*lcase(c)
char c[];
{
  int i;
   
  i = 0;
  while (c[i] != '\0') {
     if (c[i]>='A' && c[i] <= 'Z')
        c[i] = c[i] + 'a' -'A';
      i++;
       }
  return(c);
}


void
SetFont(FontName, Weight, Orientation, Size)
char *FontName, *Weight, *Orientation;
int Size;
{
   char *pattern, *instxt, SizeName[3], *LFontName, *LWeight, *LOrientation, 
        **names, *FirstName, Single[2];
   int maxnames;
   int ActCntNum;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (FontName[0] == 't' || FontName[0] == 'T')
     {
       LFontName = (char *) calloc (1, sizeof("times"));
       LFontName = strcpy(LFontName,"times");
     }
   else if (FontName[0] == 'n' || FontName[0] == 'N')
     {
       LFontName = (char *) calloc (1, sizeof("newcenturyschlbk"));
       LFontName = strcpy(LFontName,"newcenturyschlbk");
     }
   else if (FontName[0] == 'h' || FontName[0] == 'H')
     {
       LFontName = (char *) calloc (1, sizeof("helvetica"));
       LFontName = strcpy(LFontName,"helvetica");
     }
   else if (FontName[0] == 'c' || FontName[0] == 'C')
     {
       LFontName = (char *) calloc (1, sizeof("courier"));
       LFontName = strcpy(LFontName,"courier");
     }
   else if (FontName[0] == 's' || FontName[0] == 'S')
     {
       LFontName = (char *) calloc (1, sizeof("symbol"));
       LFontName = strcpy(LFontName,"symbol");
     }
   LOrientation = (char *) calloc (1, sizeof(Orientation));
   LOrientation = strcpy (LOrientation, Orientation);
   LOrientation = lcase(LOrientation);
   LWeight = (char *) calloc (1, sizeof(Weight));
   LWeight = strcpy (LWeight, Weight);
   LWeight = lcase(LWeight);
   CurrentFontName = (char *) calloc (1, 200);
   pattern = (char *) calloc (1, 200);
   CurrentFontName = strcpy (CurrentFontName, LFontName);
   CurrentFontName[0] = CurrentFontName[0]+'A'-'a';
   if (CurrentFontName[0] == 'N') 
     {
       CurrentFontName = strcpy( CurrentFontName ,"NewCenturySchlbk");
     }
   if (LOrientation[0] == 'o' && LWeight[0] == 'm')
      CurrentFontName = strcat (CurrentFontName, "-Oblique");
   else if (LOrientation[0] == 'i' && LWeight[0] == 'm')
      CurrentFontName = strcat (CurrentFontName, "-Italic");
   else if (LOrientation[0] == 'r' && LWeight[0] == 'm')
     {
       if (CurrentFontName[0] == 'T' || CurrentFontName[0] == 'N') 
         CurrentFontName = strcat (CurrentFontName, "-Roman");
     }
   else if (LOrientation[0] == 'o' && LWeight[0] == 'b')
      CurrentFontName = strcat (CurrentFontName, "-BoldOblique");
   else if (LOrientation[0] == 'i' && LWeight[0] == 'b')
      CurrentFontName = strcat (CurrentFontName, "-BoldItalic");
   else if (LOrientation[0] == 'r' && LWeight[0] == 'b')
      CurrentFontName = strcat (CurrentFontName, "-Bold");
   
   if (Output == PSOnly || Output == XAndPS)
     { 
       CurrentFontSize = (int) Size*Scale;
       fprintf (PSFile, "/%s findfont %d scalefont setfont\n", CurrentFontName,
                CurrentFontSize);
     }
   if (Output == XOnly || Output == XAndPS)
     { 
       pattern = strcpy(pattern,"-adobe-");
       strcat (pattern, LFontName); 
       strcat (pattern, "-");
       if (LWeight[0] == 'm' || LWeight[0] == 'M')
          strcat (pattern, "medium");
       else
          strcat (pattern, "bold");
       strcat (pattern, "-");
       Single[0] = Orientation[0];
       Single[1] = '\0';
       strcat (pattern, Single);
       strcat (pattern,  "-normal--");
       sprintf (SizeName, "%d", Size);
       strcat (pattern, SizeName);
       strcat (pattern, "-*");
       maxnames = 1;
       names = XListFonts(display, pattern, maxnames, &ActCntNum); 
       if (names != NULL)
         {
           FirstName = names[0];
           if (CurrentFontStruct != NULL) 
              XFreeFont (display, CurrentFontStruct);
           CurrentFontStruct=XLoadQueryFont(display, FirstName);
           CurrentFont = CurrentFontStruct->fid;
           XSetFont (display, GraphContext, CurrentFont);
           XSetFont (display, BackContext, CurrentFont);
         }
       else
         {
            printf ("Font requested not found by x-server\n");
            printf ("Requested font was: %s \n", pattern);
            FirstName = "-adobe-times-medium-r-normal--24-240-75-75-p-124-iso8859-1";
           CurrentFontStruct=XLoadQueryFont(display, FirstName);
         }
     }
   cfree (pattern);
   cfree (LFontName);
   cfree (LWeight);
   cfree (LOrientation);
   cfree(CurrentFontName);
}



void
Redraw()
{
   XSetWindowBackgroundPixmap (display, GraphWin, Backup);
   MapStartButton();
   MapStopButton();
}
  


void 
StartRedraw()

{
   XClearWindow (display, GraphWin);
   MapStartButton();
   MapStopButton();
   Redraw ();
}



void
WaitRelease()
{
   XNextEvent (display, &report);
   if (report.type == ButtonRelease)
     {
       XCloseDisplay(display);
       exit(1);
     }
}




void
CheckEvent()
{
   if (XPending(display) != 0)
     {
        XNextEvent (display, &report);
         switch (report.type)
            { 
               case (Expose)    : StartRedraw();
                                  XGetGeometry (display, GraphWin, &root, 
                                        &WinULX, &WinULY, &WinWidth,
                                        &WinHeight, &BorderWidth, &Dep); 
                                  XClearWindow (display, GraphWin);
                                  Redraw ();
                                  break;
               case ButtonPress : if ((report.xbutton.window) == StopButton) 
                          {
                            if (Output == PSOnly || Output == XAndPS)
                              { 
                                fprintf (PSFile, "stroke\ngrestore\nend\nshowpage\n"); 
                                fclose (PSFile);
                              }
                            WaitRelease();
                          }
            }
     }
}



void 
User2Win (XUser,YUser, XWin, YWin)
double XUser, YUser; 
int   *XWin, *YWin;

{
   *XWin = (XUser-UserLLX)/UserWidth*WinWidth;
   *YWin = (1.0-(YUser-UserLLY)/UserHeight)*WinHeight;
}



void 
UserPSWin (XUser,YUser, XPS, YPS)
double XUser, YUser; 
int   *XPS, *YPS;

{
   *XPS = (XUser-UserLLX)/UserWidth*PSWidth;
   *YPS = (YUser-UserLLY)/UserHeight*PSHeight;
}



void 
Win2User (XWin, YWin, XUser, YUser)
double *XUser, *YUser; 
int   XWin, YWin;

{
   *XUser = UserLLX+XWin*UserWidth/WinWidth;
   *YUser = UserLLY+(WinHeight-YWin)*UserHeight/WinHeight;
}

void
CheckPathLength(incr, X, Y)
int incr, X, Y;
{
   PathLength= PathLength+incr;
   if (PathLength>MAXPATH) 
     {
        fprintf (PSFile, "stroke\n");
        fprintf (PSFile, "%d %d M\n", X, Y);
        PathLength = 0;
     }
}




void
Draw (X1, Y1, X2, Y2)
double X1, Y1, X2, Y2;

{
   int WinX1, WinY1, WinX2, WinY2;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if ((X1<UserLLX || Y1<UserLLY || X1>UserLLX+UserWidth ||  Y1>UserLLY+UserHeight) &&
       (X2<UserLLX || Y2<UserLLY || X2>UserLLX+UserWidth ||  Y2>UserLLY+UserHeight))
     {
       if (Output == XOnly || Output == XAndPS)
         {
            User2Win (X1, Y1, &WinX1, &WinY1);
            User2Win (X2, Y2, &WinX2, &WinY2);
            LastWinX = WinX2;
            LastWinY = WinY2;
         }
       return;
     }       
   if (Output == PSOnly || Output == XAndPS)
     { 
       UserPSWin (X1, Y1, &WinX1, &WinY1);
       UserPSWin (X2, Y2, &WinX2, &WinY2);
       fprintf (PSFile, "%d %d M\n%d %d L\n", WinX1, WinY1, WinX2, WinY2);
       CheckPathLength(2, WinX2, WinY2);
     }
   if (Output == XOnly || Output == XAndPS)
     { 
       CheckEvent();
       User2Win (X1, Y1, &WinX1, &WinY1);
       User2Win (X2, Y2, &WinX2, &WinY2);

       XDrawLine (display, GraphWin, GraphContext, WinX1, WinY1, WinX2, WinY2);
       XDrawLine (display, Backup, BackContext, WinX1, WinY1, WinX2, WinY2);

       LastWinX = WinX2;
       LastWinY = WinY2;
     }
}



void 
SetPoint (X, Y)
double X, Y;
{
   int WinX, WinY;

   if (X<UserLLX || Y<UserLLY || X>UserLLX+UserWidth ||  Y>UserLLY+UserHeight)
     {
       if (Output == XOnly || Output == XAndPS)
         {
            User2Win (X, Y, &WinX, &WinY);
            LastWinX = WinX;
            LastWinY = WinY;
         }
       return;
     }       
   if (Output == PSOnly || Output == XAndPS)
     { 
       UserPSWin (X, Y, &WinX, &WinY);
       fprintf (PSFile, "%d %d M %d %d L\n", WinX, WinY, WinX, WinY);
       CheckPathLength(1, WinX, WinY);
     }
   if (Output == XOnly || Output == XAndPS)
     { 
       CheckEvent();
       User2Win (X, Y, &WinX, &WinY);
       XDrawPoint (display, GraphWin, GraphContext, WinX, WinY);
       XDrawPoint (display, Backup, BackContext, WinX, WinY);
       LastWinX = WinX;
       LastWinY = WinY;
     }
}



void 
SetDPoint (X, Y)
int X, Y;
{
   if (Output == PSOnly || Output == XAndPS)
     { 
       fprintf (PSFile, "%d %d M %d %d L\n", X, Y, X, Y);
       CheckPathLength(1, X, Y);
     }
   if (Output == XOnly || Output == XAndPS)
     { 
       CheckEvent();
       XDrawPoint (display, GraphWin, GraphContext, X, Y);
       XDrawPoint (display, Backup, BackContext, X, Y);
       LastWinX = X;
       LastWinY = Y;
     }
}




void 
DrawTo (X, Y)
double X, Y;
{
   int WinX, WinY;

   if (Output == PSOnly || Output == XAndPS)
     { 
       UserPSWin (X, Y, &WinX, &WinY);
       fprintf (PSFile, "%d %d L\n", WinX, WinY);
       CheckPathLength(1, WinX, WinY);
     }
   if (Output == XOnly || Output == XAndPS)
     { 
       CheckEvent();
       User2Win (X, Y, &WinX, &WinY);
       XDrawLine (display, GraphWin, GraphContext, LastWinX, LastWinY, 
                                                       WinX,     WinY);
       XDrawLine (display, Backup, BackContext,  LastWinX, LastWinY, 
                                                       WinX,     WinY);
       LastWinX = WinX;
       LastWinY = WinY;
     }
}

void
DrawRectangle(X1, Y1, X2, Y2)
double X1, Y1, X2, Y2;
{
   int WinX1, WinY1, WinX2, WinY2;
   unsigned int RectWidth, RectHeight;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS)
     { 
        CheckEvent();
        User2Win (X1, Y1, &WinX1, &WinY1);
        User2Win (X2, Y2, &WinX2, &WinY2);
        LastWinX = WinX1;
        LastWinY = WinY1;

        if (WinX1<WinX2)
          RectWidth = WinX2-WinX1;
        else
          {
            RectWidth = WinX1 - WinX2;
            WinX1 = WinX2;
          }
          
        if (WinY1<WinY2)
          RectHeight = WinY2-WinY1;
        else
          {
            RectHeight = WinY1 - WinY2;
            WinY1 = WinY2;
          }

        XDrawRectangle (display, GraphWin, GraphContext, WinX1, WinY1, 
                        RectWidth, RectHeight);
        XDrawRectangle (display, Backup, BackContext, WinX1, WinY1, 
                        RectWidth, RectHeight);
     }

   if (Output == PSOnly || Output == XAndPS)
     { 
       UserPSWin (X1, Y1, &WinX1, &WinY1);
       UserPSWin (X2, Y2, &WinX2, &WinY2);
       fprintf (PSFile, "%d %d M\n%d %d L\n%d %d L\n%d %d L\n%d %d L\n", WinX1, WinY1, 
                WinX2, WinY1, WinX2, WinY2, WinX1, WinY2, WinX1, WinY1);
       CheckPathLength(5, WinX1, WinY1);
     }
}


void
FillDirect(X, Y)
int X, Y;
{
  CheckEvent();
  XFillRectangle (display, GraphWin, GraphContext, X, Y, 5, 5); 
  XFillRectangle (display, Backup, BackContext, X, Y, 5, 5);
}
   



void
FillRectangle (X1, Y1, X2, Y2)
double X1, Y1, X2, Y2;
{
   int WinX1, WinY1, WinX2, WinY2;
   unsigned int RectWidth, RectHeight;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS)
     { 
       CheckEvent();
       User2Win (X1, Y1, &WinX1, &WinY1);
       User2Win (X2, Y2, &WinX2, &WinY2);
       LastWinX = WinX1;
       LastWinY = WinY1;

        if (WinX1<WinX2)
          RectWidth = WinX2-WinX1;
        else
          {
            RectWidth = WinX1 - WinX2;
            WinX1 = WinX2;
          }
          
        if (WinY1<WinY2)
          RectHeight = WinY2-WinY1;
        else
          {
            RectHeight = WinY1 - WinY2;
            WinY1 = WinY2;
          }

        XFillRectangle (display, GraphWin, GraphContext, WinX1, WinY1, 
                        RectWidth, RectHeight);
        XFillRectangle (display, Backup, BackContext, WinX1, WinY1, 
                        RectWidth, RectHeight);
      }

   if (Output == PSOnly || Output == XAndPS)
     { 
       UserPSWin (X1, Y1, &WinX1, &WinY1);
       UserPSWin (X2, Y2, &WinX2, &WinY2);
       fprintf (PSFile, "stroke newpath\n %d %d M\n%d %d L\n%d %d L\n%d %d L\n%d %d L\n", 
            WinX1, WinY1, WinX2, WinY1, WinX2, WinY2, WinX1, WinY2, WinX1, WinY1);
       fprintf (PSFile, "CLP gsave 0.00 setgray fill grestore\n");
     }
}




void FillPolygon (points, npoints)
double points[]; 
int npoints;
{
   XPoint *WinPoints;
   double UserX, UserY;
   int xx, yy;
   int i;
   
   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS) 
     {
       CheckEvent();
       WinPoints = (XPoint *) calloc (npoints, sizeof(XPoint));
     }
   if (Output == PSOnly || Output == XAndPS) 
      fprintf (PSFile, "stroke newpath\n");
   for (i=0; i<npoints; i++)
     {
       UserX = points[2*i]; 
       UserY = points[2*i+1]; 
       if (Output == XOnly || Output == XAndPS) 
         {
           User2Win (UserX, UserY, &xx, &yy); 
           WinPoints[i].x = xx;
           WinPoints[i].y = yy;
         }
       if (Output == PSOnly || Output == XAndPS) 
         {
           UserPSWin (UserX, UserY, &xx, &yy); 
           if (i==0) 
             fprintf (PSFile, "%d %d M\n", xx, yy);
           else
             fprintf (PSFile, "%d %d L\n", xx, yy);
         }
     }
   if (Output == PSOnly || Output == XAndPS) 
     {
       fprintf (PSFile, "CLP gsave 0.00 setgray fill grestore\n");
       CheckPathLength(npoints, xx, yy);
     }
   if (Output == XOnly || Output == XAndPS) 
     {
        XFillPolygon (display, GraphWin, GraphContext, WinPoints, npoints, 
                 Complex, CoordModeOrigin);
        XFillPolygon (display, Backup, BackContext, WinPoints, npoints, 
                 Complex, CoordModeOrigin);
     }
   cfree (WinPoints);
}



void
DrawCircle(X, Y, Radius)
double X, Y, Radius;
{
  double Xr, Yr, Xmax, DeltaX, YRad,
         CircPoints[10][2], UserZero;
  int PointNum =9, i, WinX, WinY, WinRad, Width, Height;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
  if (Output == XOnly || Output == XAndPS) 
    {
      CheckEvent();
      User2Win (X, Y, &WinX, &WinY);
      LastWinX = WinX;
      LastWinY = WinY;
      Xr = X+Radius;
      User2Win (Xr, Y, &WinRad, &WinY);
      WinRad = WinRad-WinX;
      WinX = WinX-WinRad;
      WinY = WinY-WinRad;
      Width = 2*WinRad;
      Height = 2*WinRad;
      XDrawArc(display, GraphWin, GraphContext, WinX, WinY, Width, Height, 0, 
               360*64-1);
      XDrawArc(display, Backup, BackContext, WinX, WinY, Width, Height, 0, 
               360*64-1);
    }
  if (Output == PSOnly || Output == XAndPS) 
    {
      UserPSWin (X, Y, &WinX, &WinY);
      Xr = X+Radius;
      UserPSWin (Xr, Y, &WinRad, &WinY);
      WinRad = WinRad - WinX;
      fprintf (PSFile, "stroke newpath %d %d %d% d %d arc\n", WinX, WinY, 
                                                    WinRad, 0, 360);
      CheckPathLength(1, WinX, WinY);
    }
}  


void 
MoveTo (X, Y)
double X, Y;
{
   int WinX, WinY;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS) 
     {
       User2Win (X, Y, &WinX, &WinY);
       LastWinX = WinX;
       LastWinY = WinY;
     }
  if (Output == PSOnly || Output == XAndPS) 
    {
      UserPSWin (X, Y, &WinX, &WinY);
      fprintf (PSFile, "%d %d M\n", WinX, WinY);
      CheckPathLength(1, WinX, WinY);
    }
}

void 
MoveRel (DX, DY)
double DX, DY;
{
  double X, Y;
  Win2User (LastWinX, LastWinY, &X, &Y);
  MoveTo (X+DX, Y+DY);
}
   
   


void 
WriteText(Text)
char *Text;
{
  int Xwin, Ywin, TextWidth;
  unsigned long returnval, CharHeight;
  Atom atom;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
  if (Output == XOnly || Output == XAndPS) 
    {
      CheckEvent();
      XGetFontProperty (CurrentFontStruct, XA_CAP_HEIGHT, &returnval); 
      CharHeight = returnval; 
      Ywin = LastWinY + CharHeight/2;
      Xwin = LastWinX;
      XDrawString (display, GraphWin, GraphContext, Xwin, Ywin, Text,
                   strlen(Text));
      XDrawString (display, Backup, GraphContext, Xwin, Ywin, Text,
                   strlen(Text));
      TextWidth = XTextWidth (CurrentFontStruct, Text, strlen(Text));
      LastWinX = LastWinX+TextWidth;
    }
  if (Output == PSOnly || Output == XAndPS) 
    fprintf (PSFile, "(%s) Lshow\n", Text);
}



void
SetDashLength (Length)
int Length;
{
  CurrentDashLength = Length;
}

void
SetLineStyle(Style)
int Style;
{
   char DashList[4];
   int PSDashLength, WinX, WinY;
   double X, Y;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS)
     { 
       if (Style == 0) 
         {
           XSetLineAttributes(display, GraphContext, CurrentLineWidth, 
                              LineSolid, CapButt, JoinMiter);
           XSetLineAttributes(display, BackContext, CurrentLineWidth, 
                              LineSolid, CapButt, JoinMiter);
           CurrentLineStyle = LineSolid;
           if (Output == XAndPS)
             fprintf (PSFile, "stroke [] 0 setdash\n");
         }
       else if (Style == 1)
         {
           DashList[0] = 2*CurrentDashLength;
           DashList[1] = CurrentDashLength;
           XSetDashes(display, GraphContext, 0, DashList, 2);
           XSetDashes(display, BackContext, 0, DashList, 2);
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           if (Output == XAndPS)
             fprintf (PSFile, "stroke \n[%d %d] 0 setdash\n", 2*PSDashLength, PSDashLength);
          } 
       else if (Style == 2)
         {
           DashList[0] = CurrentLineWidth;
           DashList[1] = CurrentDashLength;
           XSetDashes(display, GraphContext, 0, DashList, 2);
           XSetDashes(display, BackContext, 0, DashList, 2);
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           if (Output == XAndPS)
             fprintf (PSFile, "stroke \n[%d %d] 0 setdash\n", CurrentLineWidth, PSDashLength);
          } 
       else if (Style == 3)
         {
           DashList[0] = 2*CurrentDashLength;
           DashList[1] = CurrentDashLength;
           DashList[2] = CurrentLineWidth;
           DashList[3] = CurrentDashLength;
           XSetDashes(display, GraphContext, 0, DashList, 4);
           XSetDashes(display, BackContext, 0, DashList, 4);
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           if (Output == XAndPS)
             fprintf (PSFile, "stroke \n[%d %d %d %d] 0 setdash\n", 2*PSDashLength, 
                    PSDashLength, CurrentLineWidth, PSDashLength);
          } 
       if (Style != 0) 
          {
             XSetLineAttributes(display, GraphContext, CurrentLineWidth,
                                LineOnOffDash, CapButt, JoinMiter);
             XSetLineAttributes(display, BackContext, CurrentLineWidth, 
                                LineOnOffDash, CapButt, JoinMiter);
             CurrentLineStyle = LineOnOffDash;
          }
     }
   if (Output == PSOnly)
     { 
       if (Style == 0) 
         {
           CurrentLineStyle = LineSolid;
          fprintf (PSFile, "stroke [] 0 setdash\n");
         }
       else if (Style == 1)
         {
           DashList[0] = 2*CurrentDashLength;
           DashList[1] = CurrentDashLength;
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           fprintf (PSFile, "stroke \n[%d %d] 0 setdash\n", 2*PSDashLength, PSDashLength);
          } 
       else if (Style == 2)
         {
           DashList[0] = CurrentLineWidth;
           DashList[1] = CurrentDashLength;
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           fprintf (PSFile, "stroke \n[%d %d] 0 setdash\n", CurrentLineWidth, PSDashLength);
          } 
       else if (Style == 3)
         {
           DashList[0] = 2*CurrentDashLength;
           DashList[1] = CurrentDashLength;
           DashList[2] = CurrentLineWidth;
           DashList[3] = CurrentDashLength;
           PSDashLength = CurrentDashLength*PSWidth/WinWidth;
           fprintf (PSFile, "stroke \n[%d %d %d %d] 0 setdash\n", 2*PSDashLength, 
                    PSDashLength, CurrentLineWidth, PSDashLength);
          } 
       if (Style != 0) 
          {
             CurrentLineStyle = LineOnOffDash;
          }
     }
  if (Output == PSOnly || Output == XAndPS) 
  {
    Win2User(LastWinX, LastWinY, &X, &Y);
    UserPSWin(X, Y, &WinX, &WinY);
    fprintf (PSFile, "%d %d M\n", WinX, WinY);
  }
} 

void
SetLineWidth(Width)
int Width;
{
   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == XOnly || Output == XAndPS)
     { 
        XSetLineAttributes(display, GraphContext, Width, CurrentLineStyle, 
                           CapButt, JoinMiter);
        XSetLineAttributes(display, BackContext, Width, CurrentLineStyle, 
                           CapButt, JoinMiter);
        CurrentLineWidth = Width;
     }
  if (Output == PSOnly || Output == XAndPS) 
    fprintf (PSFile, "stroke %d setlinewidth\n", Width*4-3);
} 


void 
WriteFloat(Number, TotChar, DecChar)
double Number;
int TotChar, DecChar;
{
  int Xwin, Ywin, TextWidth;
  unsigned long returnval, CharHeight;
  Atom atom;
  char *format, *NewWord;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
  format = (char *) calloc (1, 20);
  NewWord = (char *) calloc (1, 64);
  strcat (format,"%");
  sprintf (NewWord, "%d.%dlf", TotChar, DecChar);
  strcat (format, NewWord);
  sprintf (NewWord, format, Number);
  if (Output == XOnly || Output == XAndPS) 
    {
      CheckEvent();
      XGetFontProperty (CurrentFontStruct, XA_CAP_HEIGHT, &returnval); 
      CharHeight = returnval; 
      Ywin = LastWinY + CharHeight/2;
      Xwin = LastWinX;
      XDrawString (display, GraphWin, GraphContext, Xwin, Ywin, NewWord,
                   strlen(NewWord));
      XDrawString (display, Backup, GraphContext, Xwin, Ywin, NewWord,
                   strlen(NewWord));
      TextWidth = XTextWidth (CurrentFontStruct, NewWord,strlen(NewWord));
      LastWinX = LastWinX+TextWidth;
    }
  if (Output == PSOnly || Output == XAndPS) 
    fprintf (PSFile, "(%s) Lshow\n", NewWord);
  cfree (format);
  cfree(NewWord);
}


void 
WriteInt(Number, CharNum)
int Number, CharNum;
{
  int Xwin, Ywin, TextWidth;
  unsigned long returnval, CharHeight;
  Atom atom;
  char *format, *NewWord;

   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
  format = (char *) calloc (1, 20);
  NewWord = (char *) calloc (1, 64);
  strcat (format, "%");
  sprintf (NewWord, "%dd", CharNum);
  strcat (format, NewWord);
  sprintf (NewWord, format, Number);
  if (Output == XOnly || Output == XAndPS) 
    {
      CheckEvent();
      XGetFontProperty (CurrentFontStruct, XA_CAP_HEIGHT, &returnval); 
      CharHeight = returnval; 
      Ywin = LastWinY + CharHeight/2;
      Xwin = LastWinX;
      XDrawString (display, GraphWin, GraphContext, Xwin, Ywin, NewWord,
                   strlen(NewWord));
      XDrawString (display, Backup, GraphContext, Xwin, Ywin, NewWord,
                   strlen(NewWord));
      TextWidth = XTextWidth (CurrentFontStruct, NewWord,strlen(NewWord));
      LastWinX = LastWinX+TextWidth;
     }
  if (Output == PSOnly || Output == XAndPS) 
       fprintf (PSFile, "(%s) Lshow\n", NewWord);
  cfree (format);
  cfree(NewWord);
}




char *
GetDirName ()
{
   char *DirName, *HomeName, LocalName;
   int i, Len;

   i = 0;
   HomeName = NULL;
   while ((environ[i] != NULL) && (strncmp("XPSDIR=", environ[i],7) != 0))
     {
       if (strncmp("HOME=", environ[i],5) == 0)
          HomeName = environ[i]+5;
       i++;
     }
   if (environ[i] == NULL)
     {
       if (HomeName != NULL)
         {
           Len = strlen(HomeName)+4;
           DirName = (char *) calloc(1, Len);
           DirName = strcpy(DirName, HomeName);
           strcat(DirName, "/XPS");
         }
       else
         DirName = ".";
     }
   else
     DirName = environ[i]+7;
   return (DirName);
}

  
#define BUFSIZE 1024

void InitPSFile(PSFileName, width,height)
unsigned int width, height;
char *PSFileName;
{
   double XSize, YSize;
   char *HdrFileName, buf[BUFSIZE], *DirName;
   int HdrFile, n, PSInt;

   WinWidth = width;
   WinHeight = height;
   PSInt = creat (PSFileName, 0644);
   DirName= GetDirName ();
   HdrFileName = (char *) calloc(1, 100);
   HdrFileName = strcpy (HdrFileName, DirName);
   strcat (HdrFileName , "/prelude");
   if ((HdrFile = open (HdrFileName, 0)) == -1)
      {
         printf ("Error: file %s not found\n", HdrFileName);   
         printf ("Use setenv XPSDIR to specify directory\n");
         printf ("containing file \'prelude\'.\n");
         exit(1);
      }
   while ((n=read(HdrFile, buf, BUFSIZE))>0)
      if (write (PSInt, buf, n) != n)
         {
            printf ("Error in copying prelude to psfile!\n");
            exit(1);
         }
   close (PSInt);
   PSFile = fopen (PSFileName, "a");
   XSize = width/9000.0;
   YSize = height/9000.0;
   if (width>height)
     {
       PSWidth = 6400;
       PSHeight = 6400*height/width;
       Scale = 6400.0/width;
       fprintf(PSFile, "%lf %lf scale\n", XSize, XSize);
     }
   else
     {
       PSWidth = 6400*width/height;
       PSHeight = 6400;
       Scale = 6400.0/height;
       fprintf(PSFile, "%lf %lf scale\n", YSize, YSize);
     }
   
   fprintf(PSFile, "90 rotate\n0 %d translate\n0 setgray\n", -PSHeight);
   fprintf(PSFile, "/Helvetica findfont 140 scalefont setfont\n");
   fprintf(PSFile, "0 0 M\n 0 %d L\n %d %d L\n %d 0 L\n CLP\n clip\n",
           PSHeight, PSWidth, PSHeight, PSWidth);
   fprintf(PSFile, "1 setlinewidth\n");
   fprintf(PSFile, "newpath\nLT0\n");
   PathLength = 0;
   cfree (HdrFileName);
}


void
PutStartButton()
{
   if (Output == XOnly || Output == XAndPS) 
     {
       StartOK = TRUE;
       StartButton = XCreateSimpleWindow (display, GraphWin, 
                                      (int) (WinWidth/60), 
                                      (int) (WinHeight/14), 40, 20, 2,
                                BlackPixel (display, screen),
                                WhitePixel (display, screen)); 
       StartContext = XCreateGC(display, StartButton, 0, NULL);
       XSetForeground (display, StartContext, BlackPixel (display, screen));
       XDrawString (display, StartButton, GraphContext, 5,12,"start",5); 
       XGrabButton (display, AnyButton, AnyModifier, StartButton, 0,
                    ButtonPressMask, GrabModeAsync, GrabModeAsync, 
                    StartButton, None);
       while (report.type != ButtonPress)
         {
            XNextEvent (display, &report);
            if (report.type == Expose)
              {
                MapStartButton();
                MapStopButton();
              }
            if ((report.xbutton.window)!=StartButton)
               report.type = Expose; 
         } 
     }
}


void
PutStopButton()
{
   if (Output == XOnly || Output == XAndPS) 
     {
       StopOK = TRUE;
       StopButton = XCreateSimpleWindow (display, GraphWin, 
                                          (int) (WinWidth/60), 
                                          (int) (WinHeight/28), 40, 20, 2, 
                                    BlackPixel (display, screen),
                                    WhitePixel (display, screen)); 
       StopContext = XCreateGC(display, StopButton, 0, NULL);
       XSetForeground (display, StopContext, BlackPixel (display, screen));
       XDrawString (display, StopButton, GraphContext, 5,12,"stop",5); 
       XGrabButton (display, AnyButton, AnyModifier, StopButton, 0,
                    ButtonPressMask|ButtonReleaseMask, GrabModeAsync, GrabModeAsync, 
                    StopButton, None);
       MapStopButton();
     }
}


void
PutContButton()
{
   if (Output == XOnly || Output == XAndPS) 
     {
       StopOK = TRUE;
       ContButton = XCreateSimpleWindow (display, GraphWin, 
                                          (int) (WinWidth/60), 
                                          (int) (WinHeight/7), 40, 20, 2, 
                                    BlackPixel (display, screen),
                                    WhitePixel (display, screen)); 
       ContContext = XCreateGC(display, ContButton, 0, NULL);
       XSetForeground (display, ContContext, BlackPixel (display, screen));
       XDrawString (display, StopButton, GraphContext, 5,12,"cont",5); 
       XGrabButton (display, AnyButton, AnyModifier, ContButton, 0,
                    ButtonPressMask|ButtonReleaseMask, GrabModeAsync, GrabModeAsync, 
                    ContButton, None);
/*       MapContButton(); */
     }
}




void
GetInfo(redmax, greenmax, bluemax, redmul, greenmul, bluemul, base)
unsigned long *redmax, *greenmax, *bluemax, *redmul, *greenmul, *bluemul, *base;
{
   mapinfo = (XStandardColormap *) malloc(sizeof(XStandardColormap)); 
        mapinfo->red_max = 7;
        mapinfo->green_max = 7;
        mapinfo->blue_max = 3;
        mapinfo->red_mult = 32;
        mapinfo->green_mult = 4;
        mapinfo->blue_mult = 1;
        mapinfo->base_pixel = 0;
  *redmax = mapinfo->red_max;
  *greenmax = mapinfo->green_max;
  *bluemax = mapinfo->blue_max;
  *redmul = mapinfo->red_mult;
  *greenmul = mapinfo->green_mult;
  *bluemul = mapinfo->blue_mult;
  *base =   mapinfo->base_pixel;
}

void
create_rgb_colormap(disp)
        Display *disp;
{
        int i, j, k, maps, base;
        XColor  *exact;
        Visual *viz;
        Window win;

/*
**      Check if the RGB_BEST_MAP resource has already been defined
*/
        disp = XOpenDisplay(NULL);
        win = RootWindow(disp, 0);
        viz = DefaultVisual(disp, 0);
/*
**      If XGetRGBColormaps returns non-zero then the colormap
**      is returned in mapinfo->colormap.
**      This can then be used as detailed in the previous example.
*/
        if (XGetRGBColormaps(disp, win, &mapinfo, &maps, XA_RGB_BEST_MAP)) {
                printf("XA_RGB_BEST_MAP is already defined \n");
                
                XCloseDisplay(disp);
                return;  
        }
/*
**      Allocate the XStandardColormap structure and create the Colormap
*/
        mapinfo = (XStandardColormap *) malloc(sizeof(XStandardColormap)); 
        mapinfo->colormap = XCreateColormap(disp, win, viz, AllocNone); 

/*
**      The following values are specific to the user's needs
*/
        mapinfo->red_max = 7;
        mapinfo->green_max = 7;
        mapinfo->blue_max = 3;
        mapinfo->red_mult = 32;
        mapinfo->green_mult = 4;
        mapinfo->blue_mult = 1;
        mapinfo->base_pixel = 0;
/*
**      These are necessary.
*/
        mapinfo->visualid = XVisualIDFromVisual(viz);
        mapinfo->killid = (XID) mapinfo->colormap;
/*
**      Now we create the RGB colour cube.
*/
        exact = (XColor *) calloc(sizeof(XColor), 256);
        base = mapinfo->base_pixel;
        for (i = 0; i < mapinfo->red_max + 1; i++) {
                for (j = 0; j < mapinfo->green_max + 1; j++) {
                        for (k = 0; k < mapinfo->blue_max + 1; k++) {
                                exact[base].blue = 65535 * k / mapinfo->blue_max;
                                exact[base].green = 65535 * j/mapinfo->green_max;
                                exact[base].red = 65535 * i / mapinfo->red_max;
                                exact[base].flags = DoRed | DoGreen | DoBlue;
                                exact[base].pixel = base++;
                        }
                }
        }
/*
**      Transfer the RGB values to the colormap
*/
        XStoreColors(disp, mapinfo->colormap, exact, 256);
/*
**      Tell the server to use this colormap as the RGB BEST MAP
*/
        XSetRGBColormaps(disp, win, mapinfo, 1, XA_RGB_BEST_MAP);
/*
**      Mark this clients resources as permanent, then return.
*/
        XSetCloseDownMode(disp, RetainPermanent);
        XCloseDisplay(disp);
        return;
}





void InitPlot(ColorName, width, height, PSFileName, OutPar)
char *ColorName, *PSFileName;
unsigned int width, height, OutPar;
{
   char *disp_name, *FirstName;
   XColor Kleur, KleurExact;
   unsigned long ColorPix;

   BackColName = (char *) calloc(1, strlen(ColorName));
   BackColName = ColorName;
   Output = OutPar;
   FramOK = StopOK = StartOK = FALSE;
   InitOK = TRUE;

   if (Output != PSOnly && Output != XOnly && Output != XAndPS) 
      error("Sorry, Output as argument to Initplot must be 0, 1 or 2");
         
   if (Output == PSOnly || Output == XAndPS) 
      InitPSFile(PSFileName, width, height);

   CurrentLineWidth = 1;
   CurrentLineStyle = LineSolid;
   CurrentDashLength = 5;
   if (Output == XOnly || Output == XAndPS)
   {
   disp_name = NULL;
/*   create_rgb_colormap(display); */
   display = XOpenDisplay (disp_name);
   if (display == NULL) 
      {
         printf ("Environment variable DISPLAY is not set. \n");
         printf ("Set this variable to hostname:display_number \n");
         printf ("where hostname is the name of the machine running \n");
         printf ("the X-server, and the display-number is usually equal \n");
         printf ("to 0.\n");
         printf ("Do not forget to add the current machine to the xhost-list!\n");
         exit(1);
      }
   screen = DefaultScreen (display);
   root = RootWindow (display, screen);
   InstallFastColors(); 
   XAllocNamedColor (display,XDefaultColormap(display,screen),ColorName, &Kleur, &KleurExact); 
   ColorPix = Kleur.pixel;
   GraphWin = XCreateSimpleWindow (display, root, 0, 0, width, height, 2, 
                                   BlackPixel (display, screen),
                                   ColorPix);
   XGetGeometry (display, GraphWin, &root, &WinULX, &WinULY, &WinWidth,
                                         &WinHeight, &BorderWidth, &Dep);
   Backup = XCreatePixmap (display, GraphWin, width, height,Dep);
   XStoreName (display, GraphWin, "Graphics Application");
   GraphContext = XCreateGC(display, GraphWin, 0, NULL);
   BackContext  = XCreateGC(display, Backup, 0, NULL);
   XSetForeground (display, GraphContext, BlackPixel (display, screen));
   XSetForeground (display, BackContext,  BlackPixel (display, screen));
   XSetBackground (display, BackContext,  ColorPix);
   XMapRaised (display, GraphWin); 
   XUngrabButton (display, AnyButton, AnyModifier, GraphWin); 
   XSelectInput (display, GraphWin, ExposureMask|ButtonPressMask|ButtonReleaseMask);
  
   SetNamedBackground(ColorName);
   XGetGeometry (display, GraphWin, &root, &WinULX, &WinULY, &WinWidth,
                                        &WinHeight, &BorderWidth, &Dep);
   Kleurtjes = XDefaultColormap(display, screen);
   CurrentFontStruct = NULL;
   FirstName = "-adobe-times-medium-r-normal--24-240-75-75-p-124-iso8859-1";
   CurrentFontStruct=XLoadQueryFont(display, FirstName);
   }
}


void
EndPlot()
{
   if (Output == PSOnly || Output == XAndPS)
     { 
       fprintf (PSFile, "stroke\ngrestore\nend\nshowpage\n"); 
       fclose (PSFile);
     }

   if (Output == PSOnly) 
     exit(1);

   if (Output == XOnly || Output == XAndPS) 
    {
      if (StopOK == FALSE) exit(1);
      while(1)
       {
         XNextEvent (display, &report);
         switch (report.type)
           { 
              case (Expose)    :  StartRedraw();
                                  XGetGeometry (display, GraphWin, &root, 
                                        &WinULX, &WinULY, &WinWidth,
                                        &WinHeight, &BorderWidth, &Dep); 
                                  XClearWindow (display, GraphWin);
                                  Redraw ();
                                  break;
              case ButtonPress : if (report.xbutton.window == StopButton)
                                   WaitRelease();
           }
      }
    }
}


void
NextPage()
{
   if (Output == PSOnly || Output == XAndPS)
     { 
       fprintf (PSFile, "stroke\ngsave\nshowpage\ngrestore\n"); 
/*       fclose (PSFile); */
     }

/*   if (Output == PSOnly) 
     exit(1); */

   if (Output == XOnly || Output == XAndPS) 
    {
      SetNamedBackground("white");
    }
}


#include <fortint.h>
