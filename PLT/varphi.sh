for dir in */*/s_* 
do 
  cd $dir
  dvips varphi1oT
  ps2eps varphi1oT.ps 
  epspdf varphi1oT.eps 
  cd -
done
#mv AIIInovm/euterpe/s_0050/varphi1oT.pdf euterpeAIII02.pdf
#mv AIIInovm/euterpe/s_0171/varphi1oT.pdf euterpeAIII04.pdf
#mv AIIInovm/euterpe/s_0357/varphi1oT.pdf euterpeAIII06.pdf
#mv AIIInovm/euterpe/s_0643/varphi1oT.pdf euterpeAIII08.pdf
#mv AIIInovm/nlambda512/s_0050/varphi1oT.pdf knososeutAIII02.pdf
#mv AIIInovm/nlambda512/s_0171/varphi1oT.pdf knososeutAIII04.pdf
#mv AIIInovm/nlambda512/s_0357/varphi1oT.pdf knososeutAIII06.pdf
#mv AIIInovm/nlambda512/s_0643/varphi1oT.pdf knososeutAIII08.pdf
#mv AIII/nlambda512/s_0050/varphi1oT.pdf knososAIII02.pdf
#mv AIII/nlambda512/s_0171/varphi1oT.pdf knososAIII04.pdf
#mv AIII/nlambda512/s_0357/varphi1oT.pdf knososAIII06.pdf
#mv AIII/nlambda512/s_0643/varphi1oT.pdf knososAIII08.pdf
mv IVnovm/euterpe/s_0047/varphi1oT.pdf euterpeIV02.pdf
mv IVnovm/euterpe/s_0172/varphi1oT.pdf euterpeIV04.pdf
mv IVnovm/euterpe/s_0359/varphi1oT.pdf euterpeIV06.pdf
mv IVnovm/euterpe/s_0641/varphi1oT.pdf euterpeIV08.pdf
mv IVnovm/nlambda256/s_0047/varphi1oT.pdf knososeutIV02.pdf
mv IVnovm/nlambda256/s_0172/varphi1oT.pdf knososeutIV04.pdf
mv IVnovm/nlambda256/s_0359/varphi1oT.pdf knososeutIV06.pdf
mv IVnovm/nlambda256/s_0641/varphi1oT.pdf knososeutIV08.pdf
mv IV/nlambda1024/s_0047/varphi1oT.pdf knososIV02.pdf
mv IV/nlambda1024/s_0172/varphi1oT.pdf knososIV04.pdf
mv IV/nlambda1024/s_0359/varphi1oT.pdf knososIV06.pdf
mv IV/nlambda1024/s_0641/varphi1oT.pdf knososIV08.pdf
