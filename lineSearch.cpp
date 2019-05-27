#include "lineSearch.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

void dcstep(double & stx,
            double & fx,
            double & dx,
            double & sty,
            double & fy,
            double & dy,
            double & stp,
            double fp,
            double dp,
            bool & brackt,
            double const stpmin,
            double const stpmax)
{
  double const zero  = 0.0,
               p66   = 0.66,
               two   = 2.0,
               three = 3.0;
  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  sgnd = dp * (dx / std::fabs(dx));
  if (fp > fx)
  {
    //    std::cout << "dcstep 1" << std::endl;
    theta = three * (fx - fp) / (stp - stx) + dx + dp;
    s = std::max(std::fabs(theta), std::max(std::fabs(dx), std::fabs(dp)));
    gamma = s * std::sqrt(std::pow(theta / s, 2) - (dx / s) * (dp / s));
    if (stp < stx) gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p / q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two) * (stp - stx);
    if (std::fabs(stpc-stx) < std::fabs(stpq-stx))
    {
      stpf = stpc;
    }
    else
    {
      stpf = stpc + (stpq - stpc)/two;
    }
    brackt = true;
  }
  else if (sgnd < zero)
  {
    //    std::cout << "dcstep 2" << std::endl;
    theta = three*(fx - fp)/(stp - stx) + dx + dp;
    s = std::max(std::fabs(theta), std::max(std::fabs(dx),std::fabs(dp)));
    gamma = s*std::sqrt(std::pow(theta/s, 2) - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (std::fabs(stpc-stp) > std::fabs(stpq-stp))
    {
      stpf = stpc;
    }
    else
    {
      stpf = stpq;
    }
    brackt = true;
  }
  else if (std::fabs(dp) < std::fabs(dx))
  {
    //    std::cout << "dcstep 3 ***************" << std::endl;
    theta = three*(fx - fp)/(stp - stx) + dx + dp;
    s = std::max(std::fabs(theta), std::max(std::fabs(dx),std::fabs(dp)));

    gamma = s*std::sqrt(std::max(zero,std::pow(theta/s, 2)-(dx/s)*(dp/s)));
    //     std::cout << "theta = " << theta << std::endl;
    //     std::cout << "gamma = " << gamma << std::endl;
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < zero && gamma != zero)
    {
      stpc = stp + r*(stx - stp);
    }
    else if (stp > stx)
    {
      stpc = stpmax;
    }
    else
    {
      stpc = stpmin;
    }
    //     std::cout << "stpc = " << stpc << std::endl;

    stpq = stp + (dp/(dp - dx))*(stx - stp);

    //     std::cout << "--------------------------------------------" << std::endl;
    //     std::cout << "stpq = " << stpq << std::endl;
    //     std::cout << "stp  = " << stp << std::endl;
    //     std::cout << "dp   = " << dp << std::endl;
    //     std::cout << "dx   = " << dx << std::endl;
    //     std::cout << "stx  = " << stx << std::endl;
    //     std::cout << "dp/(dp-dx) = " << dp/(dp-dx) << std::endl;
    //     std::cout << "--------------------------------------------" << std::endl;



    if (brackt)
    {
      if (std::fabs(stpc-stp) < std::fabs(stpq-stp))
      {
        stpf = stpc;
      }
      else
      {
        stpf = stpq;
      }
      if (stp > stx)
      {
        stpf = std::min(stp+p66*(sty-stp),stpf);
      }
      else
      {
        stpf = std::max(stp+p66*(sty-stp),stpf);
      }
    }
    else
    {
      if (std::fabs(stpc-stp) > std::fabs(stpq-stp))
      {
        stpf = stpc;
      }
      else
      {
        stpf = stpq;
      }
      stpf = std::min(stpmax,stpf);
      stpf = std::max(stpmin,stpf);
    }
  }
  else
  {
    //    std::cout << "dcstep 4" << std::endl;
    if (brackt)
    {

      theta = three*(fp - fy)/(sty - stp) + dy + dp;
      s = std::max(std::fabs(theta), std::max(std::fabs(dy),std::fabs(dp)));
      gamma = s*std::sqrt(std::pow(theta/s, 2) - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    }
    else if (stp > stx)
    {
      stpf = stpmax;
    }
    else
    {
      stpf = stpmin;
    }
  }

  if (fp > fx)
  {
    sty = stp;
    fy = fp;
    dy = dp;
  }
  else
  {
    if (sgnd < zero)
    {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  stp = stpf;
}

void dcsrch(double & f,
            double & g,
            double & stp,
            double const ftol,
            double const gtol,
            double const xtol,
            double const stpmin,
            double const stpmax,
            TaskType & task)
{
  double const zero = 0.0,
               p5   = 0.5,
               p66  = 0.66,
               xtrapl = 1.1,
               xtrapu = 4.0;
  static bool brackt;
  static int stage;
  static double finit, fx, fy, ginit,gtest, gx, gy,
         stx,sty,stmin,stmax,width,width1;
  double ftest,fm,fxm,fym, gm,gxm,gym;

  if (task == TaskType::START)
  {

//       Check the input arguments for errors.

    if (stp < stpmin)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: STP < STPMIN" << std::endl;
    }
    if (stp > stpmax)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: STP > STPMAX" << std::endl;
    }
    if (g >= zero)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: INITIAL G >= ZERO" << std::endl;
    }
    if (ftol < zero)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: FTOL < ZERO" << std::endl;
    }
    if (gtol < zero)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: GTOL < ZERO" << std::endl;
    }
    if (xtol < zero)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: XTOL < ZERO" << std::endl;
    }
    if (stpmin < zero)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: STPMIN < ZERO" << std::endl;
    }
    if (stpmax < stpmin)
    {
      task = TaskType::ERROR;
      std::cerr << "ERROR: STPMAX < STPMIN" << std::endl;
    }

//        Exit if there are errors on input.

    if (task == TaskType::ERROR)
      return;

//        Initialize local variables.

    brackt = false;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = ftol*ginit;
    width = stpmax - stpmin;
    width1 = width/p5;

//        The variables stx, fx, gx contain the values of the step,
//        function, and derivative at the best step.
//        The variables sty, fy, gy contain the value of the step,
//        function, and derivative at sty.
//        The variables stp, f, g contain the values of the step,
//        function, and derivative at stp.

    stx = zero;
    fx = finit;
    gx = ginit;
    sty = zero;
    fy = finit;
    gy = ginit;
    stmin = zero;
    stmax = stp + xtrapu*stp;
    task = TaskType::FG;

    return;

  }

//     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
//     algorithm enters the second stage.

  ftest = finit + stp*gtest;
  if (stage == 1 && f <= ftest && g >= zero)
    stage = 2;

//     Test for warnings.

  if (brackt && (stp <= stmin || stp >= stmax))
  {
    task = TaskType::WARNING;
        std::cerr << "ROUNDING ERRORS PREVENT PROGRESS\n";
  }
  if (brackt && stmax - stmin <= xtol*stmax)
  {
    task = TaskType::WARNING;
        std::cerr << "XTOL TEST SATISFIED\n";
  }
  if (stp == stpmax && f <= ftest && g <= gtest)
  {
//    task = TaskType::WARNING;
    task = TaskType::CONVERGENCE;
        std::cerr << "STP = STPMAX\n";
  }
  if (stp == stpmin && (f > ftest || g >= gtest))
  {
    task = TaskType::WARNING;
        std::cerr << "STP = STPMIN\n";
  }
//     Test for convergence.

  /*
        printf("f = %15.15g \n", f);
        printf("ftest = %15.15g\n", ftest);
        std::cout << "f " << f << std::endl;
        std::cout << "ftest " << ftest << std::endl;
        std::cout << "fabs(g) " << std::fabs(g) << std::endl;
        std::cout << "gtol*ginit " << gtol * (-ginit) << std::endl;
  */
//      if (f <= ftest && std::fabs(g) <= gtol * (-ginit))
  double eps = 1E-6;
  if (f <= ftest + eps*(std::fabs(ftest) + 1) && std::fabs(g) <= gtol * (-ginit))
    task = TaskType::CONVERGENCE;

//     Test for termination.

  if (task == TaskType::WARNING || task == TaskType::CONVERGENCE)
    return;

//     A modified function is used to predict the step during the
//     first stage if a lower function value has been obtained but
//     the decrease is not sufficient.

  if (stage == 1 && f <= fx && f > ftest)
  {
//        Define the modified function and derivative values.
    fm = f - stp*gtest;
    fxm = fx - stx*gtest;
    fym = fy - sty*gtest;
    gm = g - gtest;
    gxm = gx - gtest;
    gym = gy - gtest;
//        Call dcstep to update stx, sty, and to compute the new step.
//	 std::cout << "called dcstep 1" << std::endl;
    dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax);
//        Reset the function and derivative values for f.
    fx = fxm + stx*gtest;
    fy = fym + sty*gtest;
    gx = gxm + gtest;
    gy = gym + gtest;
  }
  else
  {

//       Call dcstep to update stx, sty, and to compute the new step.
//	std::cout << "called dcstep 2" << std::endl;
    dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax);
    //	std::cout << "stx = " << stx << std::endl;
    //	std::cout << "fx  = " << fx << std::endl;
    //	std::cout << "gx  = " << gx << std::endl;
    //	std::cout << "sty = " << sty << std::endl;
    //	std::cout << "fy  = " << fy << std::endl;
    //	std::cout << "gy  = " << gy << std::endl;
    //	std::cout << "stp = " << stp << std::endl;
    //	std::cout << "f   = " << f << std::endl;
    //	std::cout << "g   = " << g << std::endl;
    //	std::cout << "stmin = " << stmin << std::endl;
    //	std::cout << "stmax = " << stmax << std::endl;
  }
//     Decide if a bisection step is needed.

  if (brackt)
  {
    //	std::cout << "brackt true" << std::endl;
    if (std::fabs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx);
    width1 = width;
    width = std::fabs(sty-stx);
  }

//     Set the minimum and maximum steps allowed for stp.

  if (brackt)
  {
    stmin = std::min(stx,sty);
    stmax = std::max(stx,sty);
  }
  else
  {
    stmin = stp + xtrapl*(stp - stx);
    stmax = stp + xtrapu*(stp - stx);
  }
//
//     Force the step to be within the bounds stpmax and stpmin.
//
  stp = std::max(stp,stpmin);
  stp = std::min(stp,stpmax);

//     If further progress is not possible, let stp be the best
//     point obtained during the search.

  if ((brackt && (stp <= stmin || stp >= stmax))
      || (brackt && stmax-stmin <= xtol*stmax))
  {
    //	std::cout << "no further progress" << std::endl;
    stp = stx;
  }

//     Obtain another function and derivative.

  task = TaskType::FG;
}

