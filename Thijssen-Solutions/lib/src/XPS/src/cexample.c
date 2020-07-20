void main()
{
  double Rij[3][2], Size;
  int i;

  InitPlot("lightblue", 600, 600,"Tekening.ps",2);
  PutStartButton();
  PutStopButton();
  Framing (-10.0, -10.0, 10.0, 10.0);
  for (i=1; i<=18; i++)
    {
      Size = 5.+i/4.;
      DrawRectangle (-Size, -Size, Size, Size);
    } 
  SetLineStyle(1);
  SetLineWidth(1);
  Draw (-15., -5., 0., 0.);
  SetLineStyle(2);
  Draw (-15., -4., 0., 1.);
  SetLineStyle(3);
  Draw (-15., -3., 0., 2.);
  FillRectangle(-4.0, -4.0, -5.0, -5.0);
  SetFont("Times", "medium", "r", 24);
  MoveTo (0., -2.);
  WriteText ("Hallo jos");
  WriteInt (703,6);
  SetFont("Courier", "medium", "o", 18);
  MoveTo (0., 0.);
  WriteText ("Hallo jos");
  WriteFloat (703.5,6,2);
  DrawTo(0.0, 7.0);
  DrawCircle (3.0, 2.0, 5.0);
  DrawRectangle (-2.0, -3.0, 8.0, 7.0);
  SetFont ("symbol", "medium", "r", 18);
  MoveTo (0., 2.);
  WriteText ("Hallo jos");
  Rij[0][0] = 0.;
  Rij[0][1] = 0.;
  Rij[1][0] = 3.;
  Rij[1][1] = 3.;
  Rij[2][0] = 5.;
  Rij[2][1] = 0.;
  SetFillStyle("o");
  FillPolygon (Rij, 3);
  EndPlot();
}
