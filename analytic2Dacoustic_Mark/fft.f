c c********************************************************** 
c 
c 
c...          fft of real or complex array
c
c...       maximum of 8192 complex or real elements
c
c   a     =  input/output array
c   mm    =  number of elements must eqaul  2 ** mm 
c   ifr   =  1  -  forward transform
c           -1  -  inverse transform
c   irc   =  1  -  real data
c           -1  -  complex data
c
c**********************************************************
c
c
        subroutine  fft(a,mm,ifr,irc)
c
        complex  cexpt,c1,c2,c3,z,a(1)
c
c
      m = 2 ** mm
        if (irc .eq. -1)  then
        call cmpfft(a,m,ifr)
        return
        end if
c
      nn = mm - 1 
      n = m / 2 
        n1 = n + 1
        n2 = n/2 + 1
c
        if (ifr .eq. 1)  then
        call cmpfft(a,n,ifr)
        a(n1) = a(1)
        end if
          do 10 i = 1,n2
          j = n1 + 1 - i
          c1 = a(i) + conjg(a(j))
          c2 = a(i) - conjg(a(j))
          z = cexpt(i-1,m)
          c3 = cmplx(aimag(z),real(z))
          if (ifr .eq. -1)  c3 = conjg(c3)
          a(i) = 0.5 * (c1 - c2 * c3)
          a(j) = conjg(0.5 * (c1 + c2 * c3))
10        continue
        if (ifr .eq. -1)  call cmpfft(a,n,ifr)
c
        return
        end
c
c
c
c...       complex fft in place
c
	subroutine cmpfft(cx,lx,isgn)
c
	complex  cx(lx), carg, cw, ctemp, cexpt
c
	if (isgn .eq. -1)  then
		sc = 1.0 / float(lx)
			do 5 i = 1,lx
			cx(i) = sc * cx(i)
5			continue
	end if
	j = 1
	  do 30 i = 1,lx
	  if (i .gt. j)  go to 10
	  ctemp = cx(j)
	  cx(j) = cx(i)
	  cx(i) = ctemp
10	  m = lx / 2
20	  if (j .le. m)  go to 30
	  j = j - m
	  m = m / 2
	  if (m .ge. 1)  go to 20
30	  j = j + m
	l = 1
40	istep = 2 * l
	  do 50 m = 1,l
	  if (isgn .eq. 1)  then
	  cw = conjg(cexpt(m-1,l+l))
	  else
	  cw = cexpt(m-1,l+l)
	  end if
	    do 50 i = m,lx,istep
	    ctemp = cw * cx(i+l)
	    cx(i+l) = cx(i) - ctemp
50	    cx(i) = cx(i) + ctemp
	l = istep
	if (l .lt. lx)  go to 40
c
	return
	end
c
c
        complex function cexpt(k,n)
        complex z
        dimension c(2049)
c
        data  (c(i),i=    1,  128) /
     &1.0000000,.9999997,.9999988,.9999974,.9999953,.9999927,.9999894,
     & .9999856,.9999812,.9999762,.9999706,.9999644,.9999576,.9999503,
     & .9999424,.9999338,.9999247,.9999150,.9999047,.9998938,.9998823,
     & .9998703,.9998577,.9998444,.9998306,.9998162,.9998012,.9997856,
     & .9997694,.9997526,.9997353,.9997174,.9996988,.9996797,.9996600,
     & .9996397,.9996188,.9995974,.9995753,.9995527,.9995294,.9995056,
     & .9994812,.9994562,.9994306,.9994044,.9993777,.9993503,.9993224,
     & .9992939,.9992648,.9992350,.9992048,.9991739,.9991424,.9991103,
     & .9990777,.9990445,.9990107,.9989763,.9989413,.9989057,.9988695,
     & .9988328,.9987954,.9987575,.9987190,.9986799,.9986402,.9985999,
     & .9985591,.9985176,.9984756,.9984329,.9983897,.9983459,.9983016,
     & .9982566,.9982110,.9981648,.9981181,.9980708,.9980229,.9979744,
     & .9979253,.9978756,.9978253,.9977745,.9977230,.9976711,.9976184,
     & .9975652,.9975114,.9974571,.9974021,.9973466,.9972904,.9972337,
     & .9971764,.9971185,.9970601,.9970010,.9969413,.9968811,.9968203,
     & .9967589,.9966969,.9966343,.9965711,.9965074,.9964430,.9963781,
     & .9963126,.9962465,.9961798,.9961126,.9960447,.9959763,.9959072,
     & .9958376,.9957674,.9956966,.9956253,.9955533,.9954808,.9954076,
     & .9953339,.9952596/
        data  (c(i),i=  129,  256) /
     & .9951847,.9951093,.9950332,.9949566,.9948793,.9948015,.9947231,
     & .9946442,.9945646,.9944844,.9944037,.9943224,.9942405,.9941580,
     & .9940749,.9939912,.9939070,.9938222,.9937367,.9936507,.9935641,
     & .9934770,.9933892,.9933009,.9932119,.9931225,.9930323,.9929417,
     & .9928504,.9927586,.9926661,.9925731,.9924796,.9923853,.9922906,
     & .9921952,.9920993,.9920028,.9919057,.9918080,.9917098,.9916109,
     & .9915115,.9914114,.9913108,.9912097,.9911079,.9910055,.9909027,
     & .9907991,.9906950,.9905903,.9904851,.9903792,.9902728,.9901658,
     & .9900582,.9899501,.9898413,.9897320,.9896220,.9895115,.9894004,
     & .9892887,.9891765,.9890637,.9889503,.9888363,.9887217,.9886065,
     & .9884908,.9883745,.9882576,.9881401,.9880220,.9879034,.9877841,
     & .9876643,.9875439,.9874229,.9873014,.9871793,.9870566,.9869333,
     & .9868094,.9866849,.9865599,.9864343,.9863081,.9861813,.9860539,
     & .9859260,.9857975,.9856684,.9855387,.9854085,.9852777,.9851462,
     & .9850143,.9848816,.9847485,.9846148,.9844804,.9843456,.9842101,
     & .9840741,.9839374,.9838002,.9836624,.9835240,.9833851,.9832456,
     & .9831055,.9829648,.9828236,.9826817,.9825393,.9823963,.9822527,
     & .9821086,.9819639,.9818186,.9816727,.9815263,.9813792,.9812316,
     & .9810834,.9809346/
        data  (c(i),i=  257,  384) /
     & .9807853,.9806353,.9804848,.9803338,.9801821,.9800299,.9798771,
     & .9797237,.9795698,.9794152,.9792601,.9791045,.9789482,.9787914,
     & .9786339,.9784759,.9783174,.9781582,.9779985,.9778382,.9776773,
     & .9775159,.9773539,.9771913,.9770281,.9768644,.9767001,.9765352,
     & .9763697,.9762037,.9760371,.9758699,.9757021,.9755338,.9753649,
     & .9751954,.9750254,.9748547,.9746835,.9745117,.9743394,.9741665,
     & .9739929,.9738189,.9736443,.9734690,.9732932,.9731169,.9729400,
     & .9727625,.9725844,.9724057,.9722265,.9720467,.9718663,.9716854,
     & .9715039,.9713218,.9711391,.9709559,.9707721,.9705878,.9704028,
     & .9702173,.9700313,.9698446,.9696574,.9694696,.9692813,.9690923,
     & .9689028,.9687127,.9685221,.9683309,.9681391,.9679468,.9677538,
     & .9675604,.9673663,.9671717,.9669765,.9667807,.9665844,.9663875,
     & .9661900,.9659920,.9657934,.9655942,.9653944,.9651941,.9649932,
     & .9647918,.9645898,.9643872,.9641840,.9639803,.9637761,.9635712,
     & .9633658,.9631598,.9629533,.9627461,.9625385,.9623302,.9621214,
     & .9619120,.9617020,.9614916,.9612805,.9610689,.9608566,.9606438,
     & .9604305,.9602166,.9600021,.9597871,.9595715,.9593554,.9591386,
     & .9589213,.9587035,.9584851,.9582661,.9580465,.9578264,.9576057,
     & .9573845,.9571627/
        data  (c(i),i=  385,  512) /
     & .9569404,.9567174,.9564939,.9562699,.9560453,.9558201,.9555944,
     & .9553680,.9551412,.9549137,.9546857,.9544572,.9542281,.9539984,
     & .9537682,.9535374,.9533060,.9530741,.9528416,.9526086,.9523750,
     & .9521409,.9519061,.9516709,.9514350,.9511986,.9509616,.9507241,
     & .9504861,.9502474,.9500082,.9497685,.9495282,.9492873,.9490459,
     & .9488039,.9485614,.9483182,.9480746,.9478304,.9475856,.9473403,
     & .9470944,.9468479,.9466009,.9463534,.9461052,.9458566,.9456073,
     & .9453576,.9451072,.9448563,.9446048,.9443528,.9441003,.9438471,
     & .9435934,.9433392,.9430844,.9428291,.9425732,.9423168,.9420598,
     & .9418022,.9415441,.9412854,.9410262,.9407664,.9405060,.9402452,
     & .9399837,.9397218,.9394592,.9391961,.9389325,.9386683,.9384035,
     & .9381382,.9378724,.9376060,.9373390,.9370715,.9368035,.9365348,
     & .9362656,.9359959,.9357257,.9354548,.9351835,.9349116,.9346392,
     & .9343661,.9340925,.9338185,.9335437,.9332685,.9329928,.9327165,
     & .9324396,.9321622,.9318843,.9316058,.9313267,.9310471,.9307669,
     & .9304863,.9302050,.9299232,.9296409,.9293580,.9290746,.9287906,
     & .9285061,.9282210,.9279354,.9276492,.9273625,.9270753,.9267875,
     & .9264991,.9262102,.9259208,.9256308,.9253403,.9250492,.9247576,
     & .9244655,.9241728/
        data  (c(i),i=  513,  640) /
     & .9238795,.9235857,.9232914,.9229965,.9227011,.9224052,.9221087,
     & .9218116,.9215140,.9212159,.9209172,.9206180,.9203182,.9200180,
     & .9197171,.9194158,.9191139,.9188114,.9185084,.9182048,.9179008,
     & .9175962,.9172910,.9169853,.9166791,.9163723,.9160650,.9157571,
     & .9154487,.9151398,.9148303,.9145203,.9142097,.9138986,.9135870,
     & .9132749,.9129622,.9126490,.9123352,.9120209,.9117060,.9113907,
     & .9110748,.9107583,.9104413,.9101238,.9098057,.9094871,.9091680,
     & .9088483,.9085281,.9082074,.9078861,.9075643,.9072420,.9069191,
     & .9065957,.9062718,.9059473,.9056223,.9052967,.9049707,.9046441,
     & .9043170,.9039893,.9036611,.9033324,.9030031,.9026733,.9023430,
     & .9020121,.9016808,.9013488,.9010164,.9006834,.9003499,.9000159,
     & .8996813,.8993462,.8990106,.8986745,.8983378,.8980006,.8976628,
     & .8973246,.8969858,.8966464,.8963066,.8959662,.8956254,.8952839,
     & .8949420,.8945995,.8942565,.8939129,.8935689,.8932243,.8928792,
     & .8925335,.8921874,.8918407,.8914935,.8911458,.8907975,.8904487,
     & .8900994,.8897496,.8893992,.8890483,.8886970,.8883450,.8879926,
     & .8876396,.8872861,.8869321,.8865776,.8862225,.8858669,.8855109,
     & .8851542,.8847971,.8844394,.8840812,.8837225,.8833633,.8830036,
     & .8826433,.8822826/
        data  (c(i),i=  641,  768) /
     & .8819212,.8815594,.8811971,.8808343,.8804709,.8801070,.8797426,
     & .8793777,.8790122,.8786463,.8782798,.8779128,.8775453,.8771772,
     & .8768087,.8764396,.8760701,.8757000,.8753294,.8749583,.8745866,
     & .8742145,.8738418,.8734686,.8730950,.8727208,.8723460,.8719708,
     & .8715951,.8712188,.8708420,.8704648,.8700870,.8697087,.8693299,
     & .8689505,.8685707,.8681903,.8678095,.8674281,.8670462,.8666639,
     & .8662810,.8658975,.8655136,.8651292,.8647442,.8643588,.8639728,
     & .8635864,.8631994,.8628119,.8624240,.8620355,.8616465,.8612570,
     & .8608669,.8604764,.8600854,.8596938,.8593018,.8589092,.8585162,
     & .8581226,.8577286,.8573340,.8569390,.8565434,.8561473,.8557507,
     & .8553537,.8549561,.8545580,.8541594,.8537603,.8533607,.8529606,
     & .8525600,.8521589,.8517573,.8513552,.8509526,.8505495,.8501459,
     & .8497418,.8493372,.8489320,.8485264,.8481203,.8477137,.8473066,
     & .8468990,.8464909,.8460823,.8456733,.8452637,.8448536,.8444430,
     & .8440319,.8436203,.8432083,.8427957,.8423826,.8419690,.8415549,
     & .8411404,.8407254,.8403098,.8398938,.8394772,.8390602,.8386427,
     & .8382247,.8378062,.8373872,.8369677,.8365477,.8361272,.8357063,
     & .8352848,.8348629,.8344404,.8340175,.8335941,.8331702,.8327457,
     & .8323209,.8318955/
        data  (c(i),i=  769,  896) /
     & .8314696,.8310432,.8306164,.8301890,.8297612,.8293329,.8289041,
     & .8284748,.8280450,.8276148,.8271840,.8267528,.8263211,.8258888,
     & .8254561,.8250229,.8245893,.8241552,.8237205,.8232854,.8228498,
     & .8224137,.8219771,.8215401,.8211025,.8206645,.8202260,.8197870,
     & .8193475,.8189076,.8184671,.8180262,.8175848,.8171429,.8167006,
     & .8162577,.8158144,.8153706,.8149263,.8144816,.8140363,.8135906,
     & .8131444,.8126977,.8122506,.8118029,.8113548,.8109062,.8104572,
     & .8100076,.8095576,.8091071,.8086562,.8082047,.8077528,.8073004,
     & .8068476,.8063942,.8059404,.8054861,.8050313,.8045761,.8041204,
     & .8036642,.8032075,.8027504,.8022928,.8018347,.8013762,.8009171,
     & .8004577,.7999977,.7995372,.7990764,.7986150,.7981532,.7976908,
     & .7972280,.7967648,.7963011,.7958369,.7953722,.7949071,.7944415,
     & .7939755,.7935089,.7930419,.7925745,.7921066,.7916382,.7911693,
     & .7907000,.7902302,.7897599,.7892892,.7888181,.7883464,.7878743,
     & .7874017,.7869287,.7864552,.7859812,.7855068,.7850319,.7845566,
     & .7840808,.7836045,.7831278,.7826506,.7821729,.7816948,.7812163,
     & .7807372,.7802577,.7797778,.7792974,.7788165,.7783352,.7778534,
     & .7773712,.7768884,.7764053,.7759217,.7754376,.7749531,.7744681,
     & .7739827,.7734968/
        data  (c(i),i=  897, 1024) /
     & .7730104,.7725236,.7720364,.7715487,.7710605,.7705719,.7700828,
     & .7695933,.7691033,.7686129,.7681220,.7676307,.7671389,.7666467,
     & .7661540,.7656608,.7651672,.7646732,.7641788,.7636838,.7631884,
     & .7626926,.7621963,.7616996,.7612024,.7607048,.7602066,.7597082,
     & .7592092,.7587098,.7582099,.7577096,.7572088,.7567076,.7562060,
     & .7557039,.7552013,.7546984,.7541950,.7536911,.7531868,.7526820,
     & .7521768,.7516712,.7511651,.7506586,.7501516,.7496442,.7491364,
     & .7486281,.7481194,.7476102,.7471006,.7465906,.7460800,.7455691,
     & .7450578,.7445460,.7440338,.7435210,.7430080,.7424944,.7419804,
     & .7414660,.7409511,.7404358,.7399201,.7394039,.7388873,.7383703,
     & .7378528,.7373349,.7368165,.7362978,.7357786,.7352589,.7347389,
     & .7342184,.7336974,.7331761,.7326543,.7321320,.7316093,.7310863,
     & .7305627,.7300388,.7295144,.7289896,.7284644,.7279387,.7274126,
     & .7268861,.7263591,.7258317,.7253039,.7247757,.7242471,.7237180,
     & .7231885,.7226585,.7221282,.7215974,.7210662,.7205346,.7200025,
     & .7194700,.7189371,.7184038,.7178701,.7173359,.7168013,.7162663,
     & .7157308,.7151949,.7146587,.7141220,.7135848,.7130473,.7125093,
     & .7119710,.7114322,.7108930,.7103533,.7098133,.7092728,.7087319,
     & .7081906,.7076489/
        data  (c(i),i= 1025, 1152) /
     & .7071068,.7065642,.7060212,.7054778,.7049341,.7043899,.7038452,
     & .7033002,.7027547,.7022089,.7016626,.7011159,.7005688,.7000213,
     & .6994733,.6989250,.6983762,.6978270,.6972775,.6967275,.6961771,
     & .6956263,.6950751,.6945235,.6939715,.6934190,.6928661,.6923129,
     & .6917592,.6912051,.6906507,.6900958,.6895405,.6889848,.6884288,
     & .6878722,.6873153,.6867580,.6862003,.6856422,.6850836,.6845247,
     & .6839654,.6834056,.6828455,.6822850,.6817241,.6811627,.6806010,
     & .6800389,.6794763,.6789134,.6783500,.6777863,.6772221,.6766576,
     & .6760927,.6755273,.6749616,.6743955,.6738290,.6732621,.6726947,
     & .6721271,.6715589,.6709905,.6704215,.6698523,.6692826,.6687125,
     & .6681420,.6675712,.6669999,.6664283,.6658562,.6652838,.6647109,
     & .6641377,.6635641,.6629902,.6624157,.6618410,.6612658,.6606902,
     & .6601143,.6595380,.6589613,.6583841,.6578067,.6572288,.6566505,
     & .6560719,.6554928,.6549134,.6543336,.6537534,.6531729,.6525919,
     & .6520105,.6514288,.6508467,.6502641,.6496813,.6490980,.6485144,
     & .6479304,.6473460,.6467612,.6461760,.6455904,.6450045,.6444182,
     & .6438316,.6432444,.6426570,.6420692,.6414810,.6408924,.6403035,
     & .6397141,.6391245,.6385344,.6379439,.6373531,.6367618,.6361703,
     & .6355783,.6349860/
        data  (c(i),i= 1153, 1280) /
     & .6343933,.6338002,.6332067,.6326129,.6320187,.6314241,.6308292,
     & .6302339,.6296382,.6290421,.6284457,.6278490,.6272518,.6266543,
     & .6260564,.6254581,.6248595,.6242605,.6236611,.6230614,.6224613,
     & .6218608,.6212599,.6206588,.6200572,.6194553,.6188530,.6182503,
     & .6176473,.6170439,.6164402,.6158360,.6152316,.6146268,.6140215,
     & .6134160,.6128101,.6122038,.6115971,.6109902,.6103828,.6097751,
     & .6091670,.6085585,.6079498,.6073406,.6067311,.6061212,.6055110,
     & .6049004,.6042895,.6036782,.6030666,.6024546,.6018422,.6012295,
     & .6006165,.6000031,.5993893,.5987751,.5981607,.5975459,.5969307,
     & .5963151,.5956993,.5950831,.5944665,.5938495,.5932323,.5926147,
     & .5919967,.5913783,.5907597,.5901406,.5895213,.5889016,.5882815,
     & .5876611,.5870404,.5864193,.5857978,.5851761,.5845539,.5839314,
     & .5833086,.5826855,.5820619,.5814381,.5808139,.5801894,.5795645,
     & .5789393,.5783138,.5776879,.5770617,.5764351,.5758082,.5751809,
     & .5745533,.5739254,.5732971,.5726686,.5720396,.5714104,.5707807,
     & .5701508,.5695205,.5688899,.5682589,.5676277,.5669960,.5663640,
     & .5657318,.5650992,.5644662,.5638329,.5631993,.5625654,.5619311,
     & .5612965,.5606616,.5600263,.5593907,.5587547,.5581185,.5574819,
     & .5568450,.5562078/
        data  (c(i),i= 1281, 1408) /
     & .5555702,.5549323,.5542941,.5536556,.5530167,.5523775,.5517380,
     & .5510981,.5504580,.5498174,.5491766,.5485355,.5478941,.5472522,
     & .5466101,.5459677,.5453250,.5446819,.5440385,.5433948,.5427508,
     & .5421064,.5414618,.5408168,.5401714,.5395258,.5388799,.5382336,
     & .5375870,.5369402,.5362930,.5356454,.5349976,.5343494,.5337010,
     & .5330522,.5324031,.5317537,.5311040,.5304539,.5298036,.5291529,
     & .5285020,.5278507,.5271991,.5265472,.5258950,.5252425,.5245897,
     & .5239365,.5232831,.5226293,.5219753,.5213209,.5206662,.5200112,
     & .5193560,.5187004,.5180445,.5173883,.5167318,.5160750,.5154179,
     & .5147604,.5141027,.5134447,.5127864,.5121278,.5114688,.5108096,
     & .5101501,.5094903,.5088301,.5081697,.5075090,.5068479,.5061866,
     & .5055250,.5048631,.5042009,.5035384,.5028756,.5022125,.5015491,
     & .5008854,.5002214,.4995571,.4988925,.4982277,.4975625,.4968970,
     & .4962313,.4955652,.4948989,.4942323,.4935654,.4928982,.4922307,
     & .4915629,.4908948,.4902264,.4895578,.4888889,.4882196,.4875501,
     & .4868803,.4862103,.4855399,.4848692,.4841983,.4835271,.4828555,
     & .4821838,.4815117,.4808393,.4801667,.4794937,.4788205,.4781470,
     & .4774733,.4767992,.4761249,.4754502,.4747754,.4741002,.4734247,
     & .4727490,.4720730/
        data  (c(i),i= 1409, 1536) /
     & .4713967,.4707201,.4700433,.4693662,.4686888,.4680111,.4673332,
     & .4666550,.4659765,.4652977,.4646187,.4639393,.4632598,.4625799,
     & .4618998,.4612194,.4605387,.4598577,.4591765,.4584950,.4578133,
     & .4571312,.4564489,.4557664,.4550835,.4544004,.4537171,.4530334,
     & .4523496,.4516654,.4509810,.4502963,.4496113,.4489261,.4482406,
     & .4475548,.4468688,.4461825,.4454960,.4448092,.4441221,.4434348,
     & .4427472,.4420593,.4413712,.4406829,.4399942,.4393054,.4386162,
     & .4379268,.4372371,.4365472,.4358571,.4351666,.4344759,.4337850,
     & .4330938,.4324023,.4317106,.4310187,.4303265,.4296340,.4289412,
     & .4282483,.4275551,.4268616,.4261678,.4254739,.4247797,.4240852,
     & .4233904,.4226955,.4220002,.4213048,.4206091,.4199131,.4192169,
     & .4185204,.4178237,.4171267,.4164295,.4157321,.4150344,.4143365,
     & .4136383,.4129399,.4122412,.4115423,.4108431,.4101438,.4094441,
     & .4087442,.4080441,.4073438,.4066432,.4059424,.4052413,.4045400,
     & .4038384,.4031366,.4024346,.4017324,.4010299,.4003271,.3996242,
     & .3989210,.3982175,.3975139,.3968100,.3961058,.3954014,.3946968,
     & .3939920,.3932869,.3925816,.3918761,.3911704,.3904644,.3897581,
     & .3890517,.3883450,.3876381,.3869310,.3862236,.3855160,.3848082,
     & .3841002,.3833919/
        data  (c(i),i= 1537, 1664) /
     & .3826834,.3819747,.3812657,.3805566,.3798472,.3791376,.3784277,
     & .3777177,.3770074,.3762969,.3755862,.3748752,.3741640,.3734526,
     & .3727410,.3720292,.3713172,.3706049,.3698924,.3691797,.3684668,
     & .3677537,.3670403,.3663267,.3656130,.3648990,.3641848,.3634703,
     & .3627557,.3620408,.3613258,.3606105,.3598950,.3591793,.3584634,
     & .3577473,.3570309,.3563144,.3555976,.3548807,.3541635,.3534461,
     & .3527285,.3520107,.3512927,.3505745,.3498561,.3491375,.3484187,
     & .3476996,.3469804,.3462609,.3455413,.3448215,.3441014,.3433811,
     & .3426607,.3419400,.3412192,.3404981,.3397768,.3390554,.3383337,
     & .3376119,.3368898,.3361676,.3354451,.3347225,.3339996,.3332766,
     & .3325533,.3318299,.3311063,.3303824,.3296584,.3289342,.3282098,
     & .3274852,.3267604,.3260354,.3253103,.3245849,.3238593,.3231336,
     & .3224076,.3216815,.3209552,.3202287,.3195020,.3187751,.3180480,
     & .3173208,.3165933,.3158657,.3151379,.3144099,.3136817,.3129533,
     & .3122248,.3114960,.3107671,.3100380,.3093087,.3085793,.3078496,
     & .3071198,.3063897,.3056596,.3049292,.3041987,.3034679,.3027370,
     & .3020059,.3012747,.3005432,.2998116,.2990798,.2983478,.2976157,
     & .2968833,.2961509,.2954182,.2946853,.2939523,.2932191,.2924858,
     & .2917522,.2910185/
        data  (c(i),i= 1665, 1792) /
     & .2902846,.2895506,.2888164,.2880820,.2873474,.2866127,.2858778,
     & .2851427,.2844075,.2836721,.2829365,.2822008,.2814649,.2807288,
     & .2799926,.2792562,.2785197,.2777829,.2770461,.2763090,.2755718,
     & .2748344,.2740969,.2733592,.2726213,.2718833,.2711451,.2704068,
     & .2696683,.2689296,.2681908,.2674519,.2667127,.2659734,.2652340,
     & .2644944,.2637546,.2630147,.2622747,.2615345,.2607941,.2600535,
     & .2593129,.2585720,.2578311,.2570899,.2563486,.2556072,.2548656,
     & .2541239,.2533820,.2526400,.2518978,.2511554,.2504130,.2496703,
     & .2489276,.2481847,.2474416,.2466984,.2459550,.2452115,.2444679,
     & .2437241,.2429801,.2422361,.2414919,.2407475,.2400030,.2392583,
     & .2385136,.2377686,.2370236,.2362784,.2355330,.2347875,.2340419,
     & .2332962,.2325503,.2318042,.2310581,.2303118,.2295653,.2288188,
     & .2280720,.2273252,.2265782,.2258311,.2250839,.2243365,.2235890,
     & .2228414,.2220936,.2213457,.2205977,.2198495,.2191012,.2183528,
     & .2176042,.2168556,.2161068,.2153578,.2146088,.2138596,.2131103,
     & .2123608,.2116113,.2108616,.2101118,.2093619,.2086118,.2078616,
     & .2071113,.2063609,.2056104,.2048597,.2041089,.2033580,.2026070,
     & .2018559,.2011046,.2003532,.1996017,.1988501,.1980984,.1973465,
     & .1965946,.1958425/
        data  (c(i),i= 1793, 1920) /
     & .1950903,.1943380,.1935855,.1928330,.1920804,.1913276,.1905747,
     & .1898217,.1890686,.1883154,.1875621,.1868087,.1860551,.1853015,
     & .1845477,.1837938,.1830398,.1822858,.1815316,.1807773,.1800229,
     & .1792683,.1785137,.1777590,.1770042,.1762493,.1754942,.1747391,
     & .1739838,.1732285,.1724730,.1717175,.1709619,.1702061,.1694503,
     & .1686943,.1679383,.1671821,.1664259,.1656695,.1649131,.1641565,
     & .1633999,.1626432,.1618863,.1611294,.1603724,.1596153,.1588581,
     & .1581008,.1573434,.1565859,.1558284,.1550707,.1543129,.1535551,
     & .1527971,.1520391,.1512810,.1505228,.1497645,.1490061,.1482476,
     & .1474891,.1467304,.1459717,.1452129,.1444540,.1436950,.1429359,
     & .1421768,.1414175,.1406582,.1398988,.1391393,.1383797,.1376201,
     & .1368604,.1361005,.1353406,.1345807,.1338206,.1330605,.1323003,
     & .1315400,.1307796,.1300192,.1292587,.1284981,.1277374,.1269767,
     & .1262158,.1254549,.1246940,.1239329,.1231718,.1224106,.1216494,
     & .1208880,.1201266,.1193652,.1186036,.1178420,.1170803,.1163186,
     & .1155568,.1147949,.1140329,.1132709,.1125088,.1117467,.1109845,
     & .1102222,.1094598,.1086974,.1079349,.1071724,.1064098,.1056471,
     & .1048844,.1041216,.1033587,.1025958,.1018329,.1010698,.1003067,
     & .0995436,.0987804/
        data  (c(i),i= 1921, 2048) /
     & .0980171,.0972538,.0964904,.0957270,.0949635,.0941999,.0934363,
     & .0926726,.0919089,.0911451,.0903813,.0896174,.0888535,.0880895,
     & .0873255,.0865614,.0857973,.0850331,.0842688,.0835046,.0827402,
     & .0819758,.0812114,.0804469,.0796824,.0789178,.0781532,.0773885,
     & .0766238,.0758591,.0750943,.0743294,.0735645,.0727996,.0720346,
     & .0712696,.0705045,.0697394,.0689743,.0682091,.0674439,.0666786,
     & .0659133,.0651480,.0643826,.0636172,.0628517,.0620862,.0613207,
     & .0605551,.0597895,.0590239,.0582582,.0574925,.0567268,.0559610,
     & .0551952,.0544294,.0536635,.0528976,.0521317,.0513657,.0505997,
     & .0498337,.0490676,.0483016,.0475354,.0467693,.0460031,.0452369,
     & .0444707,.0437045,.0429382,.0421719,.0414056,.0406393,.0398729,
     & .0391065,.0383401,.0375736,.0368072,.0360407,.0352742,.0345077,
     & .0337411,.0329746,.0322080,.0314414,.0306748,.0299081,.0291415,
     & .0283748,.0276081,.0268414,.0260747,.0253079,.0245412,.0237744,
     & .0230076,.0222408,.0214740,.0207072,.0199404,.0191735,.0184067,
     & .0176398,.0168729,.0161061,.0153392,.0145723,.0138053,.0130384,
     & .0122715,.0115046,.0107376,.0099707,.0092037,.0084368,.0076698,
     & .0069028,.0061358,.0053689,.0046019,.0038349,.0030679,.0023009,
     & .0015339,.0007669/
        data  c(2049) /0.0/
c
        j = (8192 / n) * k
c
        if (j .le. 2048)  then
        z = cmplx(c(j+1),c(2049-j))
        go to 10
        end if
c
        if (j .le. 4096)  then
        z = cmplx(-c(4097-j),c(j-2047))
        go to 10
        end if
c
        if (j .le. 6144)  then
        z = cmplx(-c(j-4095),-c(6145-j))
        go to 10
        end if
c
        z = cmplx(c(8193-j),-c(j-6143))
c
10      cexpt = z
        return
        end

