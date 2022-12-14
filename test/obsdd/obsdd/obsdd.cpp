// obsdd.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include <stdio.h>
#include <assert.h>
#include "../../../src/rtklib.h"

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		printf("usage: obsdd obsfilename1 obsfilename2");
		return 0;
	}
	obs_t obs0 = { 0 };
	obs_t obs1 = { 0 };
	nav_t nav = { 0 };
	sta_t sta0 = { "" };
	sta_t sta1 = { "" };
	int stat, idxi = 0, idxj = 0, idxj0 = 0, i, j;
	double dt = 0.0;
	gtime_t curtime = { 0 };

	FILE *fobs = fopen("obs.txt", "w");
	FILE *fDD_G =fopen("dd_G.txt", "w");
	FILE *fDD_R = fopen("dd_R.txt", "w");
	FILE *fDD_E = fopen("dd_E.txt", "w");
	FILE *fDD_C = fopen("dd_C.txt", "w");
	FILE *fERR = fopen("dd_ERR.txt", "w");

	int ir[40] = { 0 }, iu[40] = { 0 }, ns = 0, sat[40] = { 0 }, sys[40] = { 0 }, flag[40] = { 0 }, frq[2] = { 0 };
	obsd_t obs[40+40] = { 0 };
	double resP[40][3] = { 0 };
	double resL[40][3] = { 0 };
	double curP[3] = { 0 }, curL[3] = { 0 };
	int wn = 0;
	double ws = 0.0;
	double age = 0.0;
	const int glo_frq_table[26] = { 1,
						-4,
						05,
						06,
						01,
						-4,
						05,
						06,
						-2,
						-7,
						00,
						-1,
						-2,
						-7,
						00,
						-1,
						04,
						-3,
						03,
						02,
						04,
						-3,
						03,
						02,
						 0,
						-5 };

	stat = readrnx(argv[1], 1, "", &obs0, &nav, &sta0);
	stat = readrnx(argv[2], 2, "", &obs1, &nav, &sta1);

	sortobs(&obs0);
	sortobs(&obs1);

	ns = 0;
	for (idxi = 0; idxi < obs0.n; ++idxi)
	{
		dt = timediff(obs0.data[idxi].time, curtime);
		if (curtime.time == 0.0)
		{
			curtime = obs0.data[idxi].time;
			ns = 0;
			continue;
		}
		ws = time2gpst(curtime, &wn);
		if (dt == 0.0)
		{
			for (idxj = idxj0; idxj0 < obs1.n; ++idxj)
			{
				dt = timediff(obs1.data[idxj].time, curtime);
				if (dt > 0.05)
				{
					break;
				}
				else if (dt > -0.05)
				{
					/* same epoch */
					if (obs0.data[idxi].sat == obs1.data[idxj].sat)
					{
						iu[ns] = idxi;
						ir[ns] = idxj;
						i = ns;
						j = ns + 40;
						obs[i] = obs0.data[idxi];
						obs[j] = obs1.data[idxj];
						sys[ns] = satsys(obs[i].sat, sat + ns);
						resP[ns][0] = (obs[i].P[0] == 0.0 || obs[j].P[0] == 0.0) ? (0.0) : (obs[j].P[0] - obs[i].P[0]);
						resP[ns][1] = (obs[i].P[1] == 0.0 || obs[j].P[1] == 0.0) ? (0.0) : (obs[j].P[1] - obs[i].P[1]);
						resP[ns][2] = (obs[i].P[2] == 0.0 || obs[j].P[2] == 0.0) ? (0.0) : (obs[j].P[2] - obs[i].P[2]);
						resL[ns][0] = (obs[i].L[0] == 0.0 || obs[j].L[0] == 0.0) ? (0.0) : (obs[j].L[0] - obs[i].L[0]);
						resL[ns][1] = (obs[i].L[1] == 0.0 || obs[j].L[1] == 0.0) ? (0.0) : (obs[j].L[1] - obs[i].L[1]);
						resL[ns][2] = (obs[i].L[2] == 0.0 || obs[j].L[2] == 0.0) ? (0.0) : (obs[j].L[2] - obs[i].L[2]);
						if (fobs != NULL)
						{
							fprintf(fobs, "%4i,%10.3f,%3i,%3i,%3i,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f\n", wn, ws, sat[ns], sys[ns], obs[i].rcv, obs[i].P[0], obs[i].P[1], obs[i].P[2], obs[i].L[0], obs[i].L[1], obs[i].L[2]);
							fprintf(fobs, "%4i,%10.3f,%3i,%3i,%3i,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f,%14.3f\n", wn, ws, sat[ns], sys[ns], obs[j].rcv, obs[j].P[0], obs[j].P[1], obs[j].P[2], obs[j].L[0], obs[j].L[1], obs[j].L[2]);
						}
						++ns;
					}
				}
				else
				{
					idxj0 = idxj;
				}
			}
			continue;
		}
		//
		for (i = 0; i < ns; ++i)
		{
			for (j = i + 1; j < ns; ++j)
			{
				if (sys[i] != sys[j]) continue;
				curP[0] = (resP[i][0] == 0.0 || resP[j][0] == 0.0) ? (0.0) : (resP[i][0] - resP[j][0]);
				curP[1] = (resP[i][1] == 0.0 || resP[j][1] == 0.0) ? (0.0) : (resP[i][1] - resP[j][1]);
				curP[2] = (resP[i][2] == 0.0 || resP[j][2] == 0.0) ? (0.0) : (resP[i][2] - resP[j][2]);
				curL[0] = (resL[i][0] == 0.0 || resL[j][0] == 0.0) ? (0.0) : (resL[i][0] - resL[j][0]);
				curL[1] = (resL[i][1] == 0.0 || resL[j][1] == 0.0) ? (0.0) : (resL[i][1] - resL[j][1]);
				curL[2] = (resL[i][2] == 0.0 || resL[j][2] == 0.0) ? (0.0) : (resL[i][2] - resL[j][2]);
				if (curL[0] != 0.0) curL[0] -= floor(curL[0] + 0.5);
				if (curL[1] != 0.0) curL[1] -= floor(curL[1] + 0.5);
				if (curL[2] != 0.0) curL[2] -= floor(curL[2] + 0.5);
				if (fabs(curP[0]) > 100.0 || fabs(curP[1]) > 100.0 || fabs(curP[2]) > 100.0)
				{
					fprintf(fERR, "%4i,%10.3f,%3i,%3i,%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", wn, ws, sat[i], sat[j], sys[i], curP[0], curP[1], curP[2], curL[0], curL[1], curL[2]);
				}
				else
				{
					if (sys[i] == SYS_GPS)
						fprintf(fDD_G, "%4i,%10.3f,%3i,%3i,%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", wn, ws, sat[i], sat[j], sys[i], curP[0], curP[1], curP[2], curL[0], curL[1], curL[2]);
					else if (sys[i] == SYS_GLO)
					{
						frq[0] = glo_frq_table[sat[i] - 1];
						frq[1] = glo_frq_table[sat[j] - 1];
						fprintf(fDD_R, "%4i,%10.3f,%3i,%3i,%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%3i,%3i\n", wn, ws, sat[i], sat[j], sys[i], curP[0], curP[1], curP[2], curL[0], curL[1], curL[2], frq[0] - frq[1], frq[1]);
					}
					else if (sys[i] == SYS_GAL)
						fprintf(fDD_E, "%4i,%10.3f,%3i,%3i,%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", wn, ws, sat[i], sat[j], sys[i], curP[0], curP[1], curP[2], curL[0], curL[1], curL[2]);
					else if (sys[i] == SYS_CMP)
						fprintf(fDD_C, "%4i,%10.3f,%3i,%3i,%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n", wn, ws, sat[i], sat[j], sys[i], curP[0], curP[1], curP[2], curL[0], curL[1], curL[2]);
				}
			}
		}
		//
		curtime = obs0.data[idxi].time;
		ns = 0;
	}

	free(obs0.data);
	free(obs1.data);
	free(nav.eph);
	free(nav.geph);
	free(nav.seph);

	fclose(fobs);
	fclose(fDD_G);
	fclose(fDD_R);
	fclose(fDD_E);
	fclose(fDD_C);
	fclose(fERR);
}


