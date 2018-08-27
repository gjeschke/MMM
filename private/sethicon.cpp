/*
  SETHICON - Set a new icon for the windows specified by it's handle (HWND).
  The icon is given by the HICON handle obtainedwith GETHICON.
    hicon = sethicon( hwnd, hicon, [(small=0)|big=1|all=2])

  Compiled with CPPMEX :
    cppmex sethicon user32.lib

  By Jérôme Lacaille @ Miriad Technologies (october 2002)
  http://lacaille.jerome.online.fr
  mailto:lacaille.jerome@online.fr
*/

#include "windows.h"
#include "cppmex.h"
 
/*****************************************************************************
 * The cppmex entry point.
 *****************************************************************************/
void cppMexFunction( int nlhs, mwArray mlhs[], int nrhs, const mwArray mrhs[])
{
  if ((nrhs < 2)  || !tobool(isreal(mrhs[0])) || !tobool(isreal(mrhs[1])))
    error("hicon = seticon( hwnd, hicon, [(small=0)|big=1])") ;

  int num = 0 ;
  if ((nrhs==3) && tobool(isreal(mrhs[2])))
    num = mrhs[2].Int() ;
  
  long adr = (long) mrhs[0].Double() ;
  HWND hwnd = (HWND) adr ;

  adr = (long) mrhs[1].Double() ;
  HICON hicon = (HICON) adr ;
  
  if ((num==1) || (num==2))
    hicon = (HICON)::SendMessage(hwnd, WM_SETICON, true, (LPARAM)hicon) ;
  if ((num==0) || (num==2))
    hicon = (HICON)::SendMessage(hwnd, WM_SETICON, false, (LPARAM)hicon) ;
  adr = (long) hicon ;
  mlhs[0] = mwArray(adr) ;
}

