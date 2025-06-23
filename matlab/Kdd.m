function kdd=Kdd(T)
X=1000./T;
KD_Hf=11.29 .* X - 2.275;
KD_Y= 19.47 .* X - 13.04;
KD_U= 15.32 .* X - 9.17;
KD_Th=13.02 .* X - 8.354;
Csat=ZrSaturation(T);
KD_Sm=log(13.338*Csat^(-0.622));
KD_Dy=log(2460.0*Csat^(-0.867));
KD_Yb=log(33460.*Csat^(-1.040));
KD_P= 7.646 .* X - 5.047;
kdd=exp([KD_Hf,KD_Y,KD_U,KD_Th,KD_Sm,KD_Dy,KD_Yb]); %,KD_P
end
