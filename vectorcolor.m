classdef vectorcolor
   properties (Constant)

col=[0.771624956522964,0.203215979484553,0.933482306389568;0.173635367417915,0.242037287421946,0.0177377887079281;0.114228716551283,0.390126765199177,0.649071797057105;0.00678764300280088,0.586984586668643,0.994260732010493;0.263771499061532,0.239967708290380,0.866100598142346;0.829487632677905,0.553156478555796,0.130125921805620;0.915100232580462,0.592667776113689,0.298080833676317;0.646891789519824,0.117650479206821,0.521160991438906;0.822204418506863,0.0115077508633316,0.907804505432220;0.105742220082724,0.887539524058416,0.00762165955742200;0.772520259002011,0.576157777651589,0.685264082351761;0.629456282126934,0.999989393792345,0.961102691053963;0.358460343273161,0.687621072150547,0.422044907565697;0.812167980422670,0.0264617316402420,0.708632216689409;0.806785753152915,0.217915501320807,0.0161661507763777;0.638772596721473,0.361101355849154,0.878511711562387;0.497740922766650,0.612468567891491,0.431831281949090;0.421808899684989,0.0945907079143414,0.953084645739150;0.232505679258733,0.785385022561335,0.687270816989985;0.487355846661752,0.851895850999835,0.713121492925332;0.118492597451599,0.370404267527797,0.161914378301922;0.129507224969560,0.938032094754519,0.874910707798805;0.190812979647539,0.400225730227703,0.890576318959371;0.973985805564430,0.699178348021081,0.460521143428757;0.341936759400730,0.217597213618846,0.00165367033040242;0.731067449580788,0.861279866156311,0.179350608645762;0.602064513565587,0.827812179831205,0.177028099755326;0.668257466954496,0.0891063774832089,0.124343080976103;0.470493534686605,0.895779347619325,0.268660077432229;0.527784317366683,0.318555891465713,0.844114078509106;0.588158608597368,0.674196100908368,0.723124042114830;0.0267591596623608,0.648485032291572,0.0268484851485694;0.845549675765417,0.743654245335770,0.738796537052958;0.520208920181423,0.0895878406323251,0.285154141153026;0.455818274670201,0.906020270036453,0.530269353140310;0.370329108487500,0.706765495000621,0.604297718267372;0.916705132183169,0.576251296920744,0.268732922600442;0.804354675940225,0.580848230619239,0.825556914801096;0.473877173742378,0.496288166352963,0.183831449685844;0.628526546506103,0.846973006699481,0.112525659949483;0.00408137051793345,0.960774411284413,0.272424951727263;0.943904928789095,0.278573284391079,0.0652602308723259;0.365100855305771,0.403562120182246,0.115580596182077;0.535186336404807,0.715941897359637,0.139359648899609;0.205775253869454,0.115782168027951,0.916379751142935;0.468219893149240,0.749614599545919,0.389064130870761;0.457638683532385,0.176603949871589,0.716249231648917;0.823858356607443,0.358655207879874,0.270507570373733;0.450218484646268,0.118020156300069,0.847155873933035;0.598551078361841,0.774976022604488,0.636162928636457;0.608190524537672,0.826195061099120,0.909744867555760;0.189925158675448,0.506254372307419,0.485566580499147;0.658463208818854,0.608547528860180,0.988365662644635;0.589567652387490,0.730959532737345,0.607573230300589;0.471762101738209,0.472169152963127,0.456830570648362;0.131850626390963,0.514062504427914,0.947436692904394;0.665371076928109,0.746159104539572,0.431323268110472;0.636659797928667,0.228680311999204,0.722665232981366;0.240396857042104,0.416776908927494,0.634288974226309;0.563439432042167,0.821117750668176,0.0347322996048961;0.315825033265129,0.292410966187162,0.199770533624836;0.314443194946129,0.922284174703991,0.608840581941463;0.874424592885774,0.0367976710990569,0.206006940161520;0.242975366249813,0.690559936695217,0.0257059147091947;0.569272143659410,0.446162050101890,0.170927972773452;0.366656953345114,0.611239487976973,0.290111715844084;0.0404856904087327,0.511579319034484,0.946483341481850;0.449767941163737,0.0429176795628623,0.233344206234627;0.708403788275901,0.802136671344560,0.609362124776279;0.682305062359562,0.558014257558824,0.0515911227869622;0.979921017198947,0.426326951110254,0.171523034649337;0.154655247235226,0.783318436933589,0.203847421497506;0.411519736210690,0.0778563689217342,0.478931599906270;0.821156758332884,0.908542985282074,0.625226954424668;0.818547281062126,0.685438023673637,0.904816093045064;0.0501092550256660,0.752770141010055,0.473420065909340;0.397542488888147,0.970321773242615,0.930329107568919;0.573874471386465,0.835650472162080,0.679997266259217;0.519718085328948,0.628658961580275,0.772317298408485;0.214611658461678,0.423820895734064,0.392255905181372;0.547038788716451,0.354848265031724,0.880261408219026;0.914909457211703,0.999399362148721,0.955376848547722;0.293645824844821,0.420728370056921,0.804356707119764;0.853169930977559,0.119686903119412,0.937853677951465;0.942682700138748,0.985109843314818,0.414635720911281;0.562456584135405,0.648119469889675,0.810626988374673;0.839107566085351,0.981521742258538,0.564126038846005;0.858747054507697,0.184217548231287,0.795984190345246;0.733859040518401,0.0139551822084869,0.514960203535817;0.149038680840084,0.384829165599754,0.838899970023156;0.682502149448683,0.954834232369411,0.461566812643760;0.435456179366092,0.0224548033689193,0.748298871060358;0.232479787125179,0.365360226134263,0.413214387200257;0.504194093405712,0.792057552446120,0.985665735395602;0.453115557130707,0.338969383528057,0.669959524699081;0.267604509533710,0.183723360405301,0.279211624222099;0.499964332300977,0.290177250917676,0.152579875413259;0.391385662614709,0.279297119509522,0.617017339133414;0.108939782874681,0.773848516772717,0.707066998148000;0.0902388856523021,0.299872773027414,0.749840175584137;0.0740703291035301,0.807956814839970,0.941503888118006;0.698411990479362,0.858650351620400,0.922837367442708;0.766725340597884,0.0900588926107757,0.688099473873163;0.759415311580776,0.433000021286110,0.861744883083310;0.167564421393112,0.878616533517246,0.515117865109352;0.0274352421954015,0.597624180221697,0.916679748097993;0.822655381841533,0.439832685569485,0.607275323583449;0.756258091269872,0.938248414666898,0.348387963389323;0.441195966376334,0.295106095506751,0.763951040947832;0.295469980027819,0.322267765572461,0.800074483153864;0.826651733495946,0.588896032619617,0.267284652218812;0.254931064608023,0.840566759424185,0.259321540212682;0.262733695175345,0.987732031263448,0.687148157509122;0.632736964100153,0.846928255431261,0.589847118063245;0.664080946485910,0.812661773038634,0.212397951428743;0.583029746228936,0.144254518152248,0.834673533512038;0.215017860714799,0.856142554876041,0.927629036567515;0.170265907564153,0.694598138945969,0.0189762801361323;0.0292746485601979,0.178284022887092,0.643182986811022;0.138310173415761,0.550837684329655,0.0918652684535833;0.374245824400131,0.576994160422263,0.380084070952211;0.652402835061649,0.106977770181467,0.534127963238016;0.865262668371389,0.529725912244251,0.616568054205709;0.0766414229379963,0.680494577143326,0.163564404526394;0.153659089777769,0.357355718661991,0.265044484596986;0.0922242789166398,0.141398817058500,0.147453810009356;0.0684388534763345,0.728183491567326,0.932759987337127;0.743694190789652,0.327267635343331,0.496787183500754;0.294606246848142,0.775706427503822,0.323297048152588;0.201316491195669,0.627699167945304,0.0685465616481323;0.532145286922195,0.407585808123831,0.437646164432990;0.476989686794929,0.606064149422007,0.346572372959600;0.119597204317028,0.391117404203388,0.916917143005111;0.331778283911022,0.280039554900825,0.614365968430295;0.306357856625583,0.443238631229502,0.370053942979050;0.662877989813672,0.422073398810573,0.695264765972980;0.117895419053971,0.798916682425956,0.298263151060726;0.918111037519723,0.541992798776990,0.220394920345962;0.245473668003218,0.422292840791624,0.955777641790280;0.103873257242961,0.647440707452456,0.591974710235906;0.632054087230324,0.820277730794147,0.955240851995935;0.867507701121394,0.177412447804011,0.505896488042575;0.375122457205251,0.992536156762276,0.956250548451060;0.880868500107934,0.460954061559773,0.410972845002001;0.582505578670664,0.825921885153901,0.143984476425217;0.467483813003006,0.398825252899495,0.716591882785943;0.619547377131263,0.661200553000218,0.237767327850644;0.386018388405529,0.747839947116982,0.153776726642910;0.281974530856734,0.749822109187068,0.968060775072776;0.920754716819080,0.960872670754698,0.832035015377695;0.213537125268279,0.945498111227454,0.422679461001661;0.517699115445962,0.720181761421602,0.149729618658644;0.462177162217608,0.471059558263333,0.972761554857374;0.891809394615400,0.136932632635187,0.958400603205388;0.317738123333833,0.440880762509840,0.596993897107301;0.504003693997880,0.243107804279100,0.862599702244809;0.0442662640929821,0.878146585187874,0.968139226672009;0.530882915746260,0.264664605274904,0.0995469039208802;0.453081484692742,0.296394507498748,0.737013419932354;0.494270671372603,0.271377035534227,0.606287858294789;0.610577358191280,0.756915327623808,0.943993349738325;0.710919317535448,0.685608054359800,0.758972057041176;0.837355711677154,0.723021129771069,0.801766397821380;0.667739458304444,0.628627824358881,0.400156741429779;0.384731870138558,0.241391753507502,0.0808153541741369;0.848777679792980,0.359711823322711,0.935501991820612;0.109231338424411,0.0855922747796964,0.944478726481858;0.878819424583762,0.978313333287591,0.359271676286673;0.101791606598000,0.156497817313257,0.605338018383709;0.474125603254573,0.180214040403062,0.827961131278939;0.457013943665533,0.0947389715952057,0.0351962676662402;0.219023797066548,0.983954947303213,0.657391623215351;0.222086040600375,0.654443028160903,0.0250492704043088;0.453409384579253,0.700476866533408,0.197415366989237;0.0405912841300787,0.701239119572958,0.527177854759583;0.497606457802799,0.168369321723852,0.729239371099122;0.456214468417201,0.421641002130672,0.0666341171517023;0.894678470033013,0.137604599247038,0.0722684055498000;0.114443032365886,0.711047705841404,0.362649545472327;0.793658524142241,0.645502374594721,0.588798904585788;0.390620942359674,0.473744619237287,0.237981045644375;0.217410281777132,0.275739832214688,0.994945544187353;0.716158478236135,0.808716936835129,0.555611145860300;0.171398096059648,0.532692692649645,0.247547615355969;0.722562124426943,0.730014685910713,0.165437016839451;0.301483355220672,0.597923081500905,0.603190912271892;0.962633010556508,0.0857441346646016,0.310806077918506;0.779672888838845,0.905093730489763,0.208693912242547;0.360786082360244,0.0534193733615486,0.958655133141718;0.916914877466434,0.780440605475915,0.528858089938386;0.383485503161873,0.193948727916610,0.682853251603214;0.727929790277654,0.195021986847171,0.108069399121647;0.528832281381775,0.687765718151870,0.363954881169274;0.634933903342758,0.429921635015681,0.203874569481535;0.597521544359967,0.614841602980719,0.374991674497275;0.195075707820021,0.202033821864545,0.412742115070876;0.385865843505301,0.994370075598839,0.0910230926856745;0.0467089159008134,0.409305481329933,0.408823469614839;0.270106120561086,0.0149779608159230,0.558579495646266;0.527660404902077,0.575044320152061,0.130862394314495;0.881192643721244,0.728853964817627,0.745475025213621;0.690090447687261,0.558832712951604,0.978193668401024;0.947299180847114,0.228807672911249,0.376455835251642;0.517994669490073,0.608497347825575,0.193560833061133;0.909720174529637,0.477845600566219,0.967214309620190;0.263881268723245,0.807763252453966,0.443453905596416;0.644547150970096,0.283912035879626,0.694517306775201;0.599035749961643,0.437756961368988,0.0946194599021244;0.838829048865146,0.781188779242190,0.993782248397599;0.752148840913713,0.333944475409704,0.0129505092787818;0.437512466358896,0.00913058989500137,0.196504013353541;0.386131046638966,0.00896561745026414,0.673381383160166;0.559919820357739,0.477792638090632,0.295531413935072;0.154279430355680,0.128466277461132,0.0199526660319347;0.213008589569354,0.697305248193861,0.830007917307218;0.232518290180205,0.717011422128980,0.696012887779933;0.901641691257635,0.828038654748321,0.0954915553571168;0.708052228780987,0.127152392958378,0.655480768183650;0.340826605415840,0.531251056076378,0.697205208728109;0.818793165906736,0.410276729758586,0.00329167092064020;0.893851188347549,0.600889195876178,0.344704947472164;0.780617588626156,0.790472504596030,0.863594850969367;0.818549924049973,0.226380175994546,0.731400968095302;0.0250593957126565,0.204332652712925,0.733196032579090;0.194931610475368,0.966870943000275,0.511596213080506;0.0140637262460235,0.592223537890352,0.629728649985760;0.415877500741589,0.847766484986754,0.838374368068705;0.0990084290092288,0.586170316085618,0.955847309117992;0.719197762186638,0.723921047284466,0.0223540709002636;0.998703638150778,0.599218443533811,0.817003256843155;0.650973514499598,0.0921751695148578,0.955898030461663;0.231445449588932,0.917716819214932,0.358650721988099;0.573607388677211,0.149327770096451,0.591720367840578;0.146178827814674,0.0366447925908398,0.966273860377381;0.284891018595386,0.555065404077111,0.324316331357561;0.400228097400958,0.236648038528806,0.420878520953458;0.545124414151869,0.385888775830987,0.0596621112584606;0.535194055911223,0.863735964262224,0.0381614470722691;0.163301159166638,0.458633795959286,0.0985737038900708;0.266531898014286,0.389020674162318,0.711247713939725;0.0186776326472794,0.706078782904691,0.762663203662588;0.0831579947480716,0.773172533497987,0.535191850003169;0.889605842206936,0.731749387028726,0.786650436179916;0.983548251062824,0.673149638319163,0.235760324053724;0.794289954634821,0.447012037553261,0.655113215954389;0.842199587211320,0.220903087991565,0.401516567359476;0.827683657521713,0.713013618549752,0.811577771387966;0.648756925857136,0.951614404037772,0.385455907414589;0.788247181685023,0.257374133641238,0.999850298987049;0.481646438747146,0.987195719476025,0.272364008787235;0.888985589790800,0.269503336647169,0.0227934184159022;0.492887658603488,0.526102471782751,0.835227513801047;0.740579960468909,0.958438333565065,0.272393006850393;0.529031048383610,0.424358533039598,0.668879280404516;0.198629532040597,0.600165464512451,0.696170590688232;0.176368161416673,0.777212579828001,0.899088339840030;0.0386382162673616,0.246277107055025,0.954680709304261;0.575583901772465,0.902110904504030,0.659900629047078;0.767455277677935,0.822231730643974,0.665499262429233;0.900476837956824,0.297898833351945,0.294074487177789;0.0505325608592341,0.0602644447059730,0.892499484742252;0.735393209733749,0.246619023062608,0.280592026439107;0.936770901404156,0.636092862167126,0.985332623565701;0.560115861052222,0.544529832943209,0.897322443900638;0.623012715051810,0.899398456485222,0.887575244549827;0.674862237690598,0.554849481304764,0.878788701544857;0.377788022575311,0.109716630204040,0.473485770702821;0.485642227119085,0.0643407809780867,0.581633899119483;0.715803412123794,0.923650411018820,0.217039550355672;0.898183192890998,0.449649283014034,0.711917363692362;0.364292268104546,0.428418041272029,0.584544788605248;0.219052791829439,0.810767599996887,0.112185974196059;0.782667660247255,0.230540570764669,0.477883080781676;0.401774858218534,0.405609764116425,0.204641466020115;0.788678464219070,0.341615468749730,0.267682938095621;0.234875628168866,0.0567842033007211,0.178600642859714;0.264865043164680,0.782503353965442,0.910723609407470;0.842646875692741,0.754777980134181,0.211354444580045;0.167492373941027,0.228931851316025,0.137880389813570;0.378302134058608,0.677795626999426,0.233783774885030;0.773089339525593,0.727914645909014,0.154345783350336;0.231386831899891,0.0857213232163729,0.675099103840164;0.990022977709309,0.180360721353731,0.737779495158753;0.591766485080178,0.465769919335827,0.399658172323670;0.696876112156687,0.0940773626692464,0.928539817330017;0.203055671115268,0.723396359075526,0.909968607814054;0.177087621666135,0.584568899261634,0.112099610021738;0.252882684878923,0.601051945960218,0.0618967183966571;0.235272078577302,0.958541441701207,0.891678568993048;0.479122933037640,0.174151101262959,0.598407837190304;0.705943823349398,0.680406742476591,0.838363635089045;0.970057414745866,0.569632589354736,0.335151043034308;0.696950090292747,0.431260827144043,0.211426034730042;0.271363685373030,0.777958203142642,0.727908576361571;0.716931021659451,0.709337719786380,0.547914200845408;0.927569417379861,0.354165808690560,0.495408881877805;0.974086279788169,0.508013802162745,0.199434824731547;0.0420724479538488,0.0783636512849306,0.703005107647091;0.246682655548506,0.771843015744368,0.0804575312348278;0.676268646471951,0.671989841856663,0.331976845515320;0.414378850468035,0.829107835493774,0.257089145062201;0.859434218534520,0.0940533647355800,0.509637564019902;0.420118954982075,0.00894199993731437,0.684830974052051;0.696576644004680,0.919634890954793,0.287148111772154;0.669748381702614,0.462577012171234,0.849501242344659;0.147839600661119,0.725507699325740,0.330731501611143;0.434691752719731,0.695048590613024,0.483614322303828;0.432690300367569,0.810119972810993,0.398690587515794;0.197713543123982,0.459640930359434,0.918546444301426;0.975698267469073,0.904486480010446,0.462739245267245;0.405315008309125,0.734612581942693,0.380673941963580;0.947540822699015,0.870453197540578,0.892738121436357;0.923528875485801,0.794590070879448,0.465262262557784;0.0497710271393281,0.356716707388370,0.762600875966982;0.202127826769977,0.662987069950876,0.657894382058856;0.885259747477244,0.307190382401310,0.904615677359710;0.287401533043414,0.191234908692890,0.195829571331068;0.414498044692470,0.253673641554974,0.119679259009922;0.303161130096557,0.805337273158446,0.107466173444626;0.749462934639131,0.962483606443425,0.340988940925402;0.983551828174357,0.552678589671451,0.0724531618033441;0.334263164345602,0.967825996571198,0.354671044998246;0.322313177274254,0.501703943787988,0.276803940478433;0.577681936035452,0.824837585408860,0.672131094840589;0.661101146676729,0.824667536537232,0.148971768248688;0.500332913781927,0.544130596606475,0.309307663995814;0.0360102629932049,0.00597260759940055,0.205317776002888;0.701883050939348,0.909403822213794,0.669213931781149;0.487326311903894,0.241544261755068,0.784285432936405;0.171443057793880,0.709168522999426,0.595861956161663;0.593038977160049,0.432294384213813,0.754370775752461;0.443484078237772,0.653247396889635,0.306410722161152;0.766479653347117,0.306033792876227,0.937853783002764;0.428000488858944,0.338599883987897,0.629312870703027;0.478847283480957,0.673479777901576,0.532853822537353;0.126750337733646,0.318453904429412,0.406821707318064;0.243249747858026,0.204446282833050,0.586112552108368;0.401433859165548,0.464315806130490,0.756088952717510;0.820404177337912,0.318794733259764,0.107868457442206;0.0667796222466405,0.465110543296013,0.650065333145524;0.254967289837585,0.543428957419740,0.593583660915456;0.879733199724166,0.532319250833796,0.869760228543377;0.167116984545443,0.472503173228249,0.456368624372377;0.158027507183071,0.402228712654349,0.518427462915306;0.340008904435168,0.373394362702400,0.738473219633441;0.869747167636230,0.789158628403019,0.970007389527441;0.109189726635868,0.895889664321556,0.390801642111797;0.238600361625582,0.848036793660241,0.124127514717747;0.838739515600344,0.0367164155968934,0.753286508876683;0.0331217785057532,0.674184239248582,0.296804875442746;0.0859783627974281,0.867360435327381,0.829816662674565;0.299220381177782,0.892914920623517,0.553574536756667;0.312640106829648,0.861723406551552,0.431978082873208;0.426063620003280,0.417800725663592,0.0356746279317127;0.485936694834682,0.646282271481293,0.450843324255535;0.918324131024612,0.602666684880015,0.204201117180619;0.763540055119054,0.00237310027065751,0.988481861363526;0.322794104620314,0.526860048225912,0.974322528271328;0.771696933947403,0.0366721929526350,0.00288401287005424;0.294223968349567,0.286342950189222,0.385224003593373;0.798409488847043,0.417581378063120,0.0632210100276268;0.383235665376924,0.594333737891912,0.361576147815033;0.506834198085354,0.0752400370060943,0.788261557260326;0.697564834878606,0.421637355252600,0.668426180369478;0.778780192265885,0.347833232844058,0.733139606515864;0.398657938510481,0.899565321464726,0.900962713121704;0.186675130700796,0.0133770985184927,0.615426814160621;0.617833284359247,0.551154655402261,0.241009940983526;0.853189956416626,0.688360071780602,0.875996298860962;0.620343524336981,0.232397541023085,0.404972144019678;0.694341240017909,0.949091832679241,0.671079197521697;0.675671550465186,0.771797578953578,0.706390779509270;0.288607524840775,0.278998308586543,0.942370259256444;0.491057980099600,0.410519646616264,0.862473230763139;0.798505483230020,0.409492632081116,0.0169344444362343;0.142125866445042,0.153708893304300,0.117023231268297;0.618757379988690,0.598705066841184,0.531519677420218;0.486511966356581,0.553804333023481,0.430116159761370;0.345554306242045,0.837038252083377,0.172661614964260;0.452244937333236,0.444365784460074,0.729988925468098;0.613524590062556,0.449275744377786,0.709976973617020;0.939442948799833,0.918477974129901,0.226680503343792;0.476795309019499,0.0671652037282218,0.944856851904646;0.682915976671031,0.812723790248829,0.668922252128120;0.202124955726631,0.282367823001402,0.258541107673866;0.228070522775393,0.389062497463411,0.548615598985699;0.0211049395719065,0.460794717896189,0.823494675787157;0.386600893272945,0.355468275654784,0.542294634661987;0.729486061589555,0.00545727612251490,0.741327863972311;0.764261015867801,0.780443824348761,0.864244384119165;0.670686202418718,0.496224299741534,0.521935822330146;0.0888401807958077,0.613117046174980,0.857466449175789;0.940603504649523,0.427741821708665,0.900949632406801;0.574466015806560,0.712311564824899,0.901321439269131;0.761158793964029,0.618282178014366,0.192413889965824;0.943234105438679,0.893017989365402,0.950611999316988;0.619047463694839,0.0337838463697732,0.340406695736139;0.247344051102179,0.582276867065926,0.462716355186864;0.829900421927655,0.340931086187625,0.516457096626436;0.796412899593309,0.296660471792874,0.850961304374949;0.487217906372866,0.400959518931113,0.499980158716883;0.176556054167989,0.0536732380320739,0.441920011124100;0.231924455607440,0.733182601159259,0.165759523106025;0.960539210249354,0.260638058165534,0.261740330040415;0.323662702650876,0.510519138596337,0.400705738349342;0.467562243815577,0.266790729714841,0.208868882843247;0.321268388565384,0.587316800265982,0.157647106884112;0.648272975272956,0.142889950874819,0.566501276537399;0.231123317113049,0.753565764511009,0.0891114477203482;0.0512513638093834,0.205771985801303,0.967159651986115;0.891291266008454,0.572839372263132,0.488459447958128;0.573415129154088,0.950894419935733,0.102532665876845;0.00109251880585215,0.329079339669740,0.937661931220745;0.490243249781737,0.539097154925513,0.347659847759534;0.694200770350395,0.547939619564150,0.737947106504941;0.947967103808097,0.112738254583453,0.563045802318013;0.631886170880843,0.260888076526270,0.553814581965032;0.790672639908468,0.984458789835328,0.282601057858355;0.974969429616040,0.845043144597511,0.0174848435390295;0.0907377538265451,0.209551615193795,0.314601318275730;0.783372362048606,0.866717878648546,0.0221407054807161;0.500680654508450,0.855789682729861,0.853624268250623;0.822466438427023,0.578666028929077,0.678262208447706;0.000916462051895439,0.650986538602939,0.533639813517572;0.629666842012866,0.198049205594007,0.515562874087817;0.870270150569883,0.465229356584605,0.579227766136064;0.754352903065295,0.381648063060893,0.574345970328071;0.492577207770445,0.159111444682572,0.987949795736083;0.200885132251032,0.0221243979065066,0.321315373484201;0.0251413920285640,0.468603052208940,0.0406815288664590;0.894619487543042,0.0868564711995361,0.0734099825979101;0.870868987529399,0.290583140155162,0.309542165613637;0.410555221869792,0.502317102053166,0.779014285837133;0.239697873186050,0.481222176769904,0.607144242617370;0.193222101185924,0.406870620368698,0.581212801904344;0.541026646839106,0.955146348357786,0.644551619648038;0.611215507772559,0.823471817551715,0.349467820379171;0.218954852783196,0.922164008233319,0.817401066995616;0.868571974610295,0.494016102810715,0.884338459714507;0.579391534312295,0.175634999368838,0.670446286378411;0.0786139422909350,0.884565599987344,0.992324097818407;0.143634547720818,0.336768896718628,0.953433982219476;0.561995898556246,0.558167658088921,0.977364887212717;0.980603489188012,0.669011178876040,0.871131835278661;0.383901397974172,0.0919522499400391,0.509686044223110;0.384015486463126,0.142123653736728,0.666095842864395;0.218419624835520,0.309233717460171,0.897965430319422;0.293355336639979,0.0797867021778039,0.0765564459755672;0.425392553773819,0.0617407721676284,0.299240180576160;0.323290796970121,0.544107038237489,0.120398356586532;0.0912014871755269,0.0808900681463149,0.0828994089275853;0.557413877685349,0.628740484583979,0.125891793029560;0.150521924696539,0.0921381781384300,0.863925807454378;0.0156316039782908,0.651833218176195,0.769054898024010;0.615729195341495,0.0510544115986177,0.979687952390714;0.202398896714253,0.318076239193010,0.211722100412775;0.459150581090245,0.198038313209049,0.692491435621034;0.773134432888788,0.523325084755100,0.654328169681483;0.0908189974480345,0.205184262930439,0.959354188433490;0.748371417307872,0.568919789184966,0.413552546401738;0.258844042765707,0.227427894515029,0.375434445432648;0.126486291883744,0.404098969709742,0.695763507824879;0.617718261276947,0.314959544295571,0.809080598710091;0.259498467143649,0.169172160397977,0.398676776700614;0.425471423607286,0.217009496908421,0.0637040129501248;0.937890982902197,0.255683021793886,0.887698832619640;0.286262404838413,0.377532312037757,0.700270913579268;0.279797585972381,0.734865740864959,0.465872344399364;0.263647723209654,0.345931498861314,0.654930793159702;0.0350490548559793,0.166584123627443,0.523961861226400;0.623032064872446,0.325049923825936,0.769808842392372;0.708286911798403,0.780351750676309,0.324915778911711;0.270227479690318,0.973210427506422,0.807699330716658;0.475724821631950,0.162556051159793,0.253501963990343;0.605399178226193,0.0900217353617241,0.640349141852427;0.267690237796965,0.601496125797927,0.247286947102631;0.640609297482985,0.500948077095952,0.322625129260143;0.190316001977440,0.530710205737811,0.0639949476841143;0.432744464210683,0.767478790095280,0.288143119977990;0.965699977391878,0.386843679876726,0.00231894080717299;0.635903113540234,0.0777300370104314,0.248333338045953;0.660809910078690,0.246015681502360,0.602622677082856;0.651705098053050,0.466509231062234,0.103487566970548;0.735920254595611,0.814861372274357,0.607529000108336;0.108269850632814,0.857879392194903,0.918217428602660;0.173781287290376,0.582257428590729,0.662633680291575;0.966644831813760,0.359244944700127,0.116934336402722;0.392959580213771,0.706331319407259,0.936022392958672;0.277423163503300,0.648724121652565,0.0409402608376906;0.895457468732352,0.299583063753844,0.897200928459428;0.909551063898366,0.604481305941995,0.256864368512137;0.389149125104006,0.655154992929508,0.597725413637697;0.345569669195994,0.706656938985572,0.648946162533449;0.336988967328533,0.0209941943709849,0.821704595142560;0.460465790878206,0.999529997575067,0.0691514419283276;0.596125540360714,0.383051980385398,0.00215242864302245;0.750672806745157,0.0145035879974146,0.604666740698734;0.642491079673823,0.741390763387318,0.537769456724930;0.557334073650914,0.369869301537290,0.822456979562536;0.209716686221522,0.872990066192549,0.926080050101772];

   end
end