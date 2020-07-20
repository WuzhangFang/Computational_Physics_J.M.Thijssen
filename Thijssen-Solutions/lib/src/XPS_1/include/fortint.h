
void
framing_ (Xll, Yll, Xur, Yur)
double *Xll, *Yll, *Xur, *Yur;
{
  Framing (*Xll, *Yll, *Xur, *Yur);
}

   
void
setfastcolor_ (Color)
int *Color;
{
  SetFastColor (*Color);
}



void 
setnamedcolor_ (ColorName)
char *ColorName;
{
  SetNamedColor (ColorName);
}


void
setquickcolor_(ColorPix)
int *ColorPix;
{
  SetQuickColor((unsigned long)(*ColorPix));
}

int
colorcode_(R, G, B)
int *R, *G, *B;
{
  int r;
  r = ColorCode(*R, *G, *B);
  return(r);
} 

void
setnumcolor_(ClPx)
int *ClPx;
{
  SetNumColor(ClPx);
}

void
setnamedbackground_ (ColorName)
char *ColorName;

{
   SetNamedBackground (ColorName);
}


void
setdashlength_(Length)
int *Length;
{
   SetDashLength (*Length);
}



void
setlinestyle_(Style)
int *Style;
{
   SetLineStyle(*Style);
}

setlinewidth_(Width)
int *Width;
{
   SetLineWidth(*Width);
}


void
setfillstyle_(style)
char *style;
{
  SetFillStyle(style);
}




void
setfont_(FontName, Weight, Orientation, Size)
char *FontName, *Weight, *Orientation;
int *Size;
{
   SetFont (FontName, Weight, Orientation, *Size);
}



void
draw_ (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
   Draw(*X1, *Y1, *X2, *Y2);
}



void 
setpoint_ (X, Y)
double *X, *Y;

{
   SetPoint (*X, *Y);
}


void 
drawto_ (X, Y)
double *X, *Y;
{
   DrawTo (*X, *Y);
}


void
drawrectangle_ (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  DrawRectangle (*X1, *Y1, *X2, *Y2);
}



void
fillrectangle_ (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  FillRectangle (*X1, *Y1, *X2, *Y2);
}


void 
fillpolygon_ (points, npoints)
double points[]; 
int *npoints;
{
   XPoint *WinPoints;
   double UserX, UserY;
   int xx, yy;
   int i;
   
   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == PSOnly || Output == XAndPS)
      fprintf (PSFile, "stroke newpath\n");
   if (Output == XOnly || Output == XAndPS)
     {
       CheckEvent();
       WinPoints = (XPoint *) calloc (*npoints, sizeof(XPoint));
     }
   for (i=0; i<*npoints; i++)
     {
       UserX = points[i]; 
       UserY = points[i+*npoints]; 
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
       CheckPathLength(*npoints, xx, yy);
     }
   if (Output == XOnly || Output == XAndPS)
     {
       XFillPolygon (display, GraphWin, GraphContext, WinPoints, *npoints, 
                     Complex, CoordModeOrigin);
       XFillPolygon (display, Backup, BackContext, WinPoints, *npoints, 
                     Complex, CoordModeOrigin);
     }
   cfree (WinPoints);
}





void
drawcircle_(X, Y, Radius)
double *X, *Y, *Radius;
{
   DrawCircle (*X, *Y, *Radius);
}



void
moveto_ (X, Y)
double *X, *Y;
{
  MoveTo (*X, *Y);
}



void 
writetext_(Text)
char *Text;
{
   WriteText (Text);
}


void
writefloat_(Number, TotChar, DecChar)
double *Number;
int *TotChar, *DecChar;
{
   WriteFloat(*Number, *TotChar, *DecChar);
}


void
writeint_(Number, CharNum)
int *Number,*CharNum;
{
   WriteInt(*Number, *CharNum);
}

void
putstartbutton_()
{
  PutStartButton();
}

void
putstopbutton_()
{
  PutStopButton();
}


void initplot_(ColorName, width, height, PSFileName, OutPar)
char *ColorName, *PSFileName;
unsigned int *width, *height, *OutPar;
{
   InitPlot (ColorName, *width, *height, PSFileName, *OutPar);
}

void
nextpage_ (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  NextPage ();
}

void
endplot_()
{
   EndPlot();
}


void
initbirdseye_(n1, n2, n3)
double *n1, *n2, *n3;
{
  InitBirdsEye(*n1, *n2, *n3);
}

void
d3drawline_(sx, xy, sz, ex, ey, ez)
double *sx, *xy, *sz, *ex, *ey, *ez;
{
  D3DrawLine(*sx, *xy, *sz, *ex, *ey, *ez);
}
  

void
d3setpoint_(x, y, z)
double *x, *y, *z;
{
  D3SetPoint(*x, *y, *z);
}

void
d3moveto_(x, y, z)
double *x, *y, *z;
{
  D3MoveTo(*x, *y, *z);
}


void
d3drawto_(x, y, z)
double *x, *y, *z;
{
  D3DrawTo(*x, *y, *z);
}

void
d3fillrect_(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4)
double *X1, *Y1, *Z1, *X2, *Y2, *Z2, *X3, *Y3, *Z3, *X4, *Y4, *Z4;
{
  D3FillRect(*X1, *Y1, *Z1, *X2, *Y2, *Z2, *X3, *Y3, *Z3, *X4, *Y4, *Z4);
}
/*
void
drawf(Xll, Yll, Xur, Yur, Nx, Ny)
double *Xll, *Yll, *Xur, *Yur;
int *Nx, *Ny;
{
  DrawF(*Xll, *Yll, *Xur, *Yur, *Nx, *Ny);
}

void
d3drawsurf(psi, size, maxsize)
double **psi;
int *size, *maxsize;
{
  D3DrawSurf_(psi, *size, *maxsize);
}
*/
/* Interface for HP-UX */



void
framing (Xll, Yll, Xur, Yur)
double *Xll, *Yll, *Xur, *Yur;
{
  Framing (*Xll, *Yll, *Xur, *Yur);
}

   
void
setfastcolor (Color)
int *Color;
{
  SetFastColor (*Color);
}



void 
setnamedcolor (ColorName)
char *ColorName;
{
  SetNamedColor (ColorName);
}





void
setnamedbackground (ColorName)
char *ColorName;

{
   SetNamedBackground (ColorName);
}


void
setdashlength(Length)
int *Length;
{
   SetDashLength (*Length);
}



void
setlinestyle(Style)
int *Style;
{
   SetLineStyle(*Style);
}

setlinewidth(Width)
int *Width;
{
   SetLineWidth(*Width);
}


void
setfillstyle(style)
char *style;
{
  SetFillStyle(style);
}




void
setfont(FontName, Weight, Orientation, Size)
char *FontName, *Weight, *Orientation;
int *Size;
{
   SetFont (FontName, Weight, Orientation, *Size);
}



void
draw (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
   Draw(*X1, *Y1, *X2, *Y2);
}



void 
setpoint (X, Y)
double *X, *Y;

{
   SetPoint (*X, *Y);
}


void 
drawto (X, Y)
double *X, *Y;
{
   DrawTo (*X, *Y);
}


void
drawrectangle (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  DrawRectangle (*X1, *Y1, *X2, *Y2);
}


void
fillrectangle (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  FillRectangle (*X1, *Y1, *X2, *Y2);
}


void 
fillpolygon (points, npoints)
double points[]; 
int *npoints;
{
   XPoint *WinPoints;
   double UserX, UserY;
   int xx, yy;
   int i;
   
   if (InitOK == FALSE || FramOK == FALSE) 
     error ("Sorry, you must call Initplot and Framing before anything else");
   if (Output == PSOnly || Output == XAndPS)
      fprintf (PSFile, "stroke newpath\n");
   if (Output == XOnly || Output == XAndPS)
     {
       CheckEvent();
       WinPoints = (XPoint *) calloc (*npoints, sizeof(XPoint));
     }
   for (i=0; i<*npoints; i++)
     {
       UserX = points[i]; 
       UserY = points[i+*npoints]; 
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
       CheckPathLength(*npoints, xx, yy);
     }
   if (Output == XOnly || Output == XAndPS)
     {
       XFillPolygon (display, GraphWin, GraphContext, WinPoints, *npoints, 
                     Complex, CoordModeOrigin);
       XFillPolygon (display, Backup, BackContext, WinPoints, *npoints, 
                     Complex, CoordModeOrigin);
     }
   cfree (WinPoints);
}





void
drawcircle(X, Y, Radius)
double *X, *Y, *Radius;
{
   DrawCircle (*X, *Y, *Radius);
}



void
moveto (X, Y)
double *X, *Y;
{
  MoveTo (*X, *Y);
}



void 
writetext(Text)
char *Text;
{
   WriteText (Text);
}


void
writefloat(Number, TotChar, DecChar)
double *Number;
int *TotChar, *DecChar;
{
   WriteFloat(*Number, *TotChar, *DecChar);
}


void
writeint(Number, CharNum)
int *Number,*CharNum;
{
   WriteInt(*Number, *CharNum);
}

void
putstartbutton()
{
  PutStartButton();
}

void
putstopbutton()
{
  PutStopButton();
}


void initplot(ColorName, width, height, PSFileName, OutPar)
char *ColorName, *PSFileName;
unsigned int *width, *height, *OutPar;
{
   InitPlot (ColorName, *width, *height, PSFileName, *OutPar);
}


void
nextpage (X1, Y1, X2, Y2)
double *X1, *Y1, *X2, *Y2;
{
  NextPage ();
}

void
endplot()
{
   EndPlot();
}
