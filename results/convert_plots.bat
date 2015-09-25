@ECHO OFF
SETLOCAL ENABLEDELAYEDEXPANSION

FOR /D %%G in (*) DO (
	CD %%G
	CALL ps2pdf perpontok.ps %%G.pdf
	CD ..
)
PAUSE

