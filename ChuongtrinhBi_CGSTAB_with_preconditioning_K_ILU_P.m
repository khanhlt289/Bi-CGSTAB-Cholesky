% Chuong trinh thu nghiem giai he phuong trinh dai so tuyen tinh Ax = b
% bang phuong phap Bi_CGSTAB co tien dieu kien phai K 
% Tien dieu kien phai K la 1 phan tich ILU_P, p = 13.
% Ngay thu nghiem: 9-2-2021
% Ngay sua chua: 22-2-2021 (lan 1)
% Nguoi lap trinh va thu nghiem: Luu Truong Khanh
% Muc dich: Phuc vu lam luan van tot nghiep
% Ket qua: Thu nghiem da thanh cong! Nghiem xap xi da giai duoc.
% Nguoi kiem tra: Luu Truong Khanh. Ngay kiem tra: 10-2-2021.
% Phat trien phan mem nang cao: Luu Truong Khanh
% Ket qua phat trien va sua chua (lan 1): Thanh cong!
% Chay thu nghiem voi chuong trinh san UBoot 12x12m. Ket qua: Thanh cong!
% Chay thu nghiem voi chuong trinh tuong 2 lop. Ket qua: ...
% =====================================================================
% Chon A, b
% clear all; clc;
% Lap he phuong trinh dai so tuyen tinh Ax = b
clear all; clc;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Nhap kich thuoc ma tran do cung, so nut, so chuyen vi nut.
format long % Dinh dang so
nel=840;   % So luong phan tu cua he
nnel=2;     % So nut cua 1 phan tu
ndof=6;     % So chuyen vi nut dau phan tu
nnode=441;   % So nut cua he
sdof=nnode*ndof;    % So chuyen vi nut cua he
lx = 0.6; % 
lz = 0.6; % 
L = 0.6; % 
% Kich thuoc dam 
b = 0.3; % Chieu rong dam la 30 cm.
h = 0.6; % Chieu cao dam la 60 cm.
E = 2.65*10^6;   % Mo dun dan hoi vat lieu T/m2
Iz =  b*h^3/12;
Iy = h*b^3/12;
Ix = Iz*0.75;
G = E/2; % Mo dun dan hoi truot vat lieu T/m2       
Area = 0.18;    % Dien tich tiet dien dam
n = 400; % So luong o san la 400 o (dang don gian)

for i=1:21
    x(i) = (i - 1)* 0.6;
    z(i) = 0;
    y(i) = 0;    
end
for i=22:42
    x(i) = (i - 22)* 0.6;
    z(i) = 0.6;
    y(i) = 0;    
end
for i=43:63
    x(i) = (i - 43)* 0.6;
    z(i) = 1.2;
    y(i) = 0;    
end
for i=64:84
    x(i) = (i - 64)* 0.6;
    z(i) = 1.8;
    y(i) = 0;    
end
for i=85:105
    x(i) = (i - 85)* 0.6;
    z(i) = 2.4;
    y(i) = 0;    
end
for i=106:126
    x(i) = (i - 106)* 0.6;
    z(i) = 3.0;
    y(i) = 0;    
end

for i=127:147
    x(i) = (i - 127)* 0.6;
    z(i) = 3.6;
    y(i) = 0;    
end

for i=148:168
    x(i) = (i - 148)* 0.6;
    z(i) = 4.2;
    y(i) = 0;    
end

for i=169:189
    x(i) = (i - 169)* 0.6;
    z(i) = 4.8;
    y(i) = 0;    
end

for i=190:210
    x(i) = (i - 190)* 0.6;
    z(i) = 5.4;
    y(i) = 0;    
end

for i=211:231
    x(i) = (i - 211)* 0.6;
    z(i) = 6.0;
    y(i) = 0;    
end

for i=232:252
    x(i) = (i - 232)* 0.6;
    z(i) = 6.6;
    y(i) = 0;    
end

for i=253:273
    x(i) = (i - 253)* 0.6;
    z(i) = 7.2;
    y(i) = 0;    
end

for i=274:294
    x(i) = (i - 274)* 0.6;
    z(i) = 7.8;
    y(i) = 0;    
end

for i=295:315
    x(i) = (i - 295)* 0.6;
    z(i) = 8.4;
    y(i) = 0;    
end

for i=316:336
    x(i) = (i - 316)* 0.6;
    z(i) = 9.0;
    y(i) = 0;    
end

for i=337:357
    x(i) = (i - 337)* 0.6;
    z(i) = 9.6;
    y(i) = 0;    
end

for i=358:378
    x(i) = (i - 358)* 0.6;
    z(i) = 10.2;
    y(i) = 0;    
end

for i=379:399
    x(i) = (i - 379)* 0.6;
    z(i) = 10.8;
    y(i) = 0;    
end

for i=400:420
    x(i) = (i - 400)* 0.6;
    z(i) = 11.4;
    y(i) = 0;    
end

for i=421:441
    x(i) = (i - 421)* 0.6;
    z(i) = 12.0;
    y(i) = 0;    
end

% Danh dau nhan cua cac doan thanh ngan
node(1,1)=1;node(1,2)=2;    
node(2,1)=2;node(2,2)=3;    
node(3,1)=3;node(3,2)=4;   
node(4,1)=4;node(4,2)=5;
node(5,1)=5;node(5,2)=6;
node(6,1)=6;node(6,2)=7;
node(7,1)=7;node(7,2)=8;
node(8,1)=8;node(8,2)=9;
node(9,1)=9;node(9,2)=10;
node(10,1)=10;node(10,2)=11;
node(11,1)=11;node(11,2)=12;
node(12,1)=12;node(12,2)=13;
node(13,1)=13;node(13,2)=14;   
node(14,1)=14;node(14,2)=15;    
node(15,1)=15;node(15,2)=16;   
node(16,1)=16;node(16,2)=17;
node(17,1)=17;node(17,2)=18;
node(18,1)=18;node(18,2)=19;
node(19,1)=19;node(19,2)=20;
node(20,1)=20;node(20,2)=21;

node(21,1)=1;node(21,2)=22;
node(22,1)=2;node(22,2)=23;
node(23,1)=3;node(23,2)=24;
node(24,1)=4;node(24,2)=25;
node(25,1)=5;node(25,2)=26;
node(26,1)=6;node(26,2)=27;
node(27,1)=7;node(27,2)=28;
node(28,1)=8;node(28,2)=29;
node(29,1)=9;node(29,2)=30;
node(30,1)=10;node(30,2)=31;
node(31,1)=11;node(31,2)=32;
node(32,1)=12;node(32,2)=33;
node(33,1)=13;node(33,2)=34;
node(34,1)=14;node(34,2)=35;
node(35,1)=15;node(35,2)=36;
node(36,1)=16;node(36,2)=37;
node(37,1)=17;node(37,2)=38;
node(38,1)=18;node(38,2)=39;
node(39,1)=19;node(39,2)=40;
node(40,1)=20;node(40,2)=41;
node(41,1)=21;node(41,2)=42;

node(42,1)=22;node(42,2)=23;   
node(43,1)=23;node(43,2)=24;    
node(44,1)=24;node(44,2)=25;   
node(45,1)=25;node(45,2)=26;
node(46,1)=26;node(46,2)=27;
node(47,1)=27;node(47,2)=28;
node(48,1)=28;node(48,2)=29;
node(49,1)=29;node(49,2)=30;
node(50,1)=30;node(50,2)=31;
node(51,1)=31;node(51,2)=32;
node(52,1)=32;node(52,2)=33;
node(53,1)=33;node(53,2)=34;
node(54,1)=34;node(54,2)=35;   
node(55,1)=35;node(55,2)=36;   
node(56,1)=36;node(56,2)=37;   
node(57,1)=37;node(57,2)=38;
node(58,1)=38;node(58,2)=39;
node(59,1)=39;node(59,2)=40;
node(60,1)=40;node(60,2)=41;
node(61,1)=41;node(61,2)=42;

node(62,1)=22;node(62,2)=43;
node(63,1)=23;node(63,2)=44;
node(64,1)=24;node(64,2)=45;
node(65,1)=25;node(65,2)=46;
node(66,1)=26;node(66,2)=47;
node(67,1)=27;node(67,2)=48;
node(68,1)=28;node(68,2)=49;
node(69,1)=29;node(69,2)=50;
node(70,1)=30;node(70,2)=51;
node(71,1)=31;node(71,2)=52;
node(72,1)=32;node(72,2)=53;
node(73,1)=33;node(73,2)=54;
node(74,1)=34;node(74,2)=55;
node(75,1)=35;node(75,2)=56;
node(76,1)=36;node(76,2)=57;
node(77,1)=37;node(77,2)=58;
node(78,1)=38;node(78,2)=59;
node(79,1)=39;node(79,2)=60;
node(80,1)=40;node(80,2)=61;
node(81,1)=41;node(81,2)=62;
node(82,1)=42;node(82,2)=63;

node(83,1)=43;node(83,2)=44;    
node(84,1)=44;node(84,2)=45;    
node(85,1)=45;node(85,2)=46;    
node(86,1)=46;node(86,2)=47;
node(87,1)=47;node(87,2)=48;
node(88,1)=48;node(88,2)=49;
node(89,1)=49;node(89,2)=50;
node(90,1)=50;node(90,2)=51;
node(91,1)=51;node(91,2)=52;
node(92,1)=52;node(92,2)=53;
node(93,1)=53;node(93,2)=54;
node(94,1)=54;node(94,2)=55;
node(95,1)=55;node(95,2)=56;    
node(96,1)=56;node(96,2)=57;   
node(97,1)=57;node(97,2)=58;   
node(98,1)=58;node(98,2)=59;
node(99,1)=59;node(99,2)=60;
node(100,1)=60;node(100,2)=61;
node(101,1)=61;node(101,2)=62;
node(102,1)=62;node(102,2)=63;

node(103,1) = 43; node(103,2) = 64;
node(104,1) = 44; node(104,2) = 65;
node(105,1) = 45; node(105,2) = 66;
node(106,1) = 46; node(106,2) = 67;
node(107,1) = 47; node(107,2) = 68;
node(108,1) = 48; node(108,2) = 69;
node(109,1) = 49; node(109,2) = 70;
node(110,1) = 50; node(110,2) = 71;
node(111,1) = 51; node(111,2) = 72;
node(112,1) = 52; node(112,2) = 73;
node(113,1) = 53; node(113,2) = 74;
node(114,1) = 54; node(114,2) = 75;
node(115,1) = 55; node(115,2) = 76;
node(116,1) = 56; node(116,2) = 77;
node(117,1) = 57; node(117,2) = 78;
node(118,1) = 58; node(118,2) = 79;
node(119,1) = 59; node(119,2) = 80;
node(120,1) = 60; node(120,2) = 81;
node(121,1) = 61; node(121,2) = 82;
node(122,1) = 62; node(122,2) = 83;
node(123,1) = 63; node(123,2) = 84;

node(124,1) = 64; node(124,2) = 65;
node(125,1) = 65; node(125,2) = 66;
node(126,1) = 66; node(126,2) = 67;
node(127,1) = 67; node(127,2) = 68;
node(128,1) = 68; node(128,2) = 69;
node(129,1) = 69; node(129,2) = 70;
node(130,1) = 70; node(130,2) = 71;
node(131,1) = 71; node(131,2) = 72;
node(132,1) = 72; node(132,2) = 73;
node(133,1) = 73; node(133,2) = 74;
node(134,1) = 74; node(134,2) = 75;
node(135,1) = 75; node(135,2) = 76;
node(136,1) = 76; node(136,2) = 77;
node(137,1) = 77; node(137,2) = 78;
node(138,1) = 78; node(138,2) = 79;
node(139,1) = 79; node(139,2) = 80;
node(140,1) = 80; node(140,2) = 81;
node(141,1) = 81; node(141,2) = 82;
node(142,1) = 82; node(142,2) = 83;
node(143,1) = 83; node(143,2) = 84;

node(144,1) = 64; node(144,2) = 85;
node(145,1) = 65; node(145,2) = 86;
node(146,1) = 66; node(146,2) = 87;
node(147,1) = 67; node(147,2) = 88;
node(148,1) = 68; node(148,2) = 89;
node(149,1) = 69; node(149,2) = 90;
node(150,1) = 70; node(150,2) = 91;
node(151,1) = 71; node(151,2) = 92;
node(152,1) = 72; node(152,2) = 93;
node(153,1) = 73; node(153,2) = 94;
node(154,1) = 74; node(154,2) = 95;
node(155,1) = 75; node(155,2) = 96;
node(156,1) = 76; node(156,2) = 97;
node(157,1) = 77; node(157,2) = 98;
node(158,1) = 78; node(158,2) = 99;
node(159,1) = 79; node(159,2) = 100;
node(160,1) = 80; node(160,2) = 101;
node(161,1) = 81; node(161,2) = 102;
node(162,1) = 82; node(162,2) = 103;
node(163,1) = 83; node(163,2) = 104;
node(164,1) = 84; node(164,2) = 105;

node(165,1) = 85; node(165,2) = 86;
node(166,1) = 86; node(166,2) = 87;
node(167,1) = 87; node(167,2) = 88;
node(168,1) = 88; node(168,2) = 89;
node(169,1) = 89; node(169,2) = 90;
node(170,1) = 90; node(170,2) = 91;
node(171,1) = 91; node(171,2) = 92;
node(172,1) = 92; node(172,2) = 93;
node(173,1) = 93; node(173,2) = 94;
node(174,1) = 94; node(174,2) = 95;
node(175,1) = 95; node(175,2) = 96;
node(176,1) = 96; node(176,2) = 97;
node(177,1) = 97; node(177,2) = 98;
node(178,1) = 98; node(178,2) = 99;
node(179,1) = 99; node(179,2) = 100;
node(180,1) = 100; node(180,2) = 101;
node(181,1) = 101; node(181,2) = 102;
node(182,1) = 102; node(182,2) = 103;
node(183,1) = 103; node(183,2) = 104;
node(184,1) = 104; node(184,2) = 105;

node(185,1) = 85; node(185,2) = 106;
node(186,1) = 86; node(186,2) = 107;
node(187,1) = 87; node(187,2) = 108;
node(188,1) = 88; node(188,2) = 109;
node(189,1) = 89; node(189,2) = 110;
node(190,1) = 90; node(190,2) = 111;
node(191,1) = 91; node(191,2) = 112;
node(192,1) = 92; node(192,2) = 113;
node(193,1) = 93; node(193,2) = 114;
node(194,1) = 94; node(194,2) = 115;
node(195,1) = 95; node(195,2) = 116;
node(196,1) = 96; node(196,2) = 117;
node(197,1) = 97; node(197,2) = 118;
node(198,1) = 98; node(198,2) = 119;
node(199,1) = 99; node(199,2) = 120;
node(200,1) = 100; node(200,2) = 121;
node(201,1) = 101; node(201,2) = 122;
node(202,1) = 102; node(202,2) = 123;
node(203,1) = 103; node(203,2) = 124;
node(204,1) = 104; node(204,2) = 125;
node(205,1) = 105; node(205,2) = 126;

node(206,1) = 106; node(206,2) = 107;
node(207,1) = 107; node(207,2) = 108;
node(208,1) = 108; node(208,2) = 109;
node(209,1) = 109; node(209,2) = 110;
node(210,1) = 110; node(210,2) = 111;
node(211,1) = 111; node(211,2) = 112;
node(212,1) = 112; node(212,2) = 113;
node(213,1) = 113; node(213,2) = 114;
node(214,1) = 114; node(214,2) = 115;
node(215,1) = 115; node(215,2) = 116;
node(216,1) = 116; node(216,2) = 117;
node(217,1) = 117; node(217,2) = 118;
node(218,1) = 118; node(218,2) = 119;
node(219,1) = 119; node(219,2) = 120;
node(220,1) = 120; node(220,2) = 121;
node(221,1) = 121; node(221,2) = 122;
node(222,1) = 122; node(222,2) = 123;
node(223,1) = 123; node(223,2) = 124;
node(224,1) = 124; node(224,2) = 125;
node(225,1) = 125; node(225,2) = 126;

node(226,1) = 106; node(226,2) = 127;
node(227,1) = 107; node(227,2) = 128;
node(228,1) = 108; node(228,2) = 129;
node(229,1) = 109; node(229,2) = 130;
node(230,1) = 110; node(230,2) = 131;
node(231,1) = 111; node(231,2) = 132;
node(232,1) = 112; node(232,2) = 133;
node(233,1) = 113; node(233,2) = 134;
node(234,1) = 114; node(234,2) = 135;
node(235,1) = 115; node(235,2) = 136;
node(236,1) = 116; node(236,2) = 137;
node(237,1) = 117; node(237,2) = 138;
node(238,1) = 118; node(238,2) = 139;
node(239,1) = 119; node(239,2) = 140;
node(240,1) = 120; node(240,2) = 141;
node(241,1) = 121; node(241,2) = 142;
node(242,1) = 122; node(242,2) = 143;
node(243,1) = 123; node(243,2) = 144;
node(244,1) = 124; node(244,2) = 145;
node(245,1) = 125; node(245,2) = 146;
node(246,1) = 126; node(246,2) = 147;

node(247,1) = 127; node(247,2) = 128;
node(248,1) = 128; node(248,2) = 129;
node(249,1) = 129; node(249,2) = 130;
node(250,1) = 130; node(250,2) = 131;
node(251,1) = 131; node(251,2) = 132;
node(252,1) = 132; node(252,2) = 133;
node(253,1) = 133; node(253,2) = 134;
node(254,1) = 134; node(254,2) = 135;
node(255,1) = 135; node(255,2) = 136;
node(256,1) = 136; node(256,2) = 137;
node(257,1) = 137; node(257,2) = 138;
node(258,1) = 138; node(258,2) = 139;
node(259,1) = 139; node(259,2) = 140;
node(260,1) = 140; node(260,2) = 141;
node(261,1) = 141; node(261,2) = 142;
node(262,1) = 142; node(262,2) = 143;
node(263,1) = 143; node(263,2) = 144;
node(264,1) = 144; node(264,2) = 145;
node(265,1) = 145; node(265,2) = 146;
node(266,1) = 146; node(266,2) = 147;

node(267,1) = 127; node(267,2) = 148;
node(268,1) = 128; node(268,2) = 149;
node(269,1) = 129; node(269,2) = 150;
node(270,1) = 130; node(270,2) = 151;
node(271,1) = 131; node(271,2) = 152;
node(272,1) = 132; node(272,2) = 153;
node(273,1) = 133; node(273,2) = 154;
node(274,1) = 134; node(274,2) = 155;
node(275,1) = 135; node(275,2) = 156;
node(276,1) = 136; node(276,2) = 157;
node(277,1) = 137; node(277,2) = 158;
node(278,1) = 138; node(278,2) = 159;
node(279,1) = 139; node(279,2) = 160;
node(280,1) = 140; node(280,2) = 161;
node(281,1) = 141; node(281,2) = 162;
node(282,1) = 142; node(282,2) = 163;
node(283,1) = 143; node(283,2) = 164;
node(284,1) = 144; node(284,2) = 165;
node(285,1) = 145; node(285,2) = 166;
node(286,1) = 146; node(286,2) = 167;
node(287,1) = 147; node(287,2) = 168;

node(288,1) = 148; node(288,2) = 149;
node(289,1) = 149; node(289,2) = 150;
node(290,1) = 150; node(290,2) = 151;
node(291,1) = 151; node(291,2) = 152;
node(292,1) = 152; node(292,2) = 153;
node(293,1) = 153; node(293,2) = 154;
node(294,1) = 154; node(294,2) = 155;
node(295,1) = 155; node(295,2) = 156;
node(296,1) = 156; node(296,2) = 157;
node(297,1) = 157; node(297,2) = 158;
node(298,1) = 158; node(298,2) = 159;
node(299,1) = 159; node(299,2) = 160;
node(300,1) = 160; node(300,2) = 161;
node(301,1) = 161; node(301,2) = 162;
node(302,1) = 162; node(302,2) = 163;
node(303,1) = 163; node(303,2) = 164;
node(304,1) = 164; node(304,2) = 165;
node(305,1) = 165; node(305,2) = 166;
node(306,1) = 166; node(306,2) = 167;
node(307,1) = 167; node(307,2) = 168;

node(308,1) = 148; node(308,2) = 169;
node(309,1) = 149; node(309,2) = 170;
node(310,1) = 150; node(310,2) = 171;
node(311,1) = 151; node(311,2) = 172;
node(312,1) = 152; node(312,2) = 173;
node(313,1) = 153; node(313,2) = 174;
node(314,1) = 154; node(314,2) = 175;
node(315,1) = 155; node(315,2) = 176;
node(316,1) = 156; node(316,2) = 177;
node(317,1) = 157; node(317,2) = 178;
node(318,1) = 158; node(318,2) = 179;
node(319,1) = 159; node(319,2) = 180;
node(320,1) = 160; node(320,2) = 181;
node(321,1) = 161; node(321,2) = 182;
node(322,1) = 162; node(322,2) = 183;
node(323,1) = 163; node(323,2) = 184;
node(324,1) = 164; node(324,2) = 185;
node(325,1) = 165; node(325,2) = 186;
node(326,1) = 166; node(326,2) = 187;
node(327,1) = 167; node(327,2) = 188;
node(328,1) = 168; node(328,2) = 189;

node(329,1) = 169; node(329,2) = 170;
node(330,1) = 170; node(330,2) = 171;
node(331,1) = 171; node(331,2) = 172;
node(332,1) = 172; node(332,2) = 173;
node(333,1) = 173; node(333,2) = 174;
node(334,1) = 174; node(334,2) = 175;
node(335,1) = 175; node(335,2) = 176;
node(336,1) = 176; node(336,2) = 177;
node(337,1) = 177; node(337,2) = 178;
node(338,1) = 178; node(338,2) = 179;
node(339,1) = 179; node(339,2) = 180;
node(340,1) = 180; node(340,2) = 181;
node(341,1) = 181; node(341,2) = 182;
node(342,1) = 182; node(342,2) = 183;
node(343,1) = 183; node(343,2) = 184;
node(344,1) = 184; node(344,2) = 185;
node(345,1) = 185; node(345,2) = 186;
node(346,1) = 186; node(346,2) = 187;
node(347,1) = 187; node(347,2) = 188;
node(348,1) = 188; node(348,2) = 189;

node(349,1) = 169; node(349,2) = 190;
node(350,1) = 170; node(350,2) = 191;
node(351,1) = 171; node(351,2) = 192;
node(352,1) = 172; node(352,2) = 193;
node(353,1) = 173; node(353,2) = 194;
node(354,1) = 174; node(354,2) = 195;
node(355,1) = 175; node(355,2) = 196;
node(356,1) = 176; node(356,2) = 197;
node(357,1) = 177; node(357,2) = 198;
node(358,1) = 178; node(358,2) = 199;
node(359,1) = 179; node(359,2) = 200;
node(360,1) = 180; node(360,2) = 201;
node(361,1) = 181; node(361,2) = 202;
node(362,1) = 182; node(362,2) = 203;
node(363,1) = 183; node(363,2) = 204;
node(364,1) = 184; node(364,2) = 205;
node(365,1) = 185; node(365,2) = 206;
node(366,1) = 186; node(366,2) = 207;
node(367,1) = 187; node(367,2) = 208;
node(368,1) = 188; node(368,2) = 209;
node(369,1) = 189; node(369,2) = 210;

node(370,1) = 190; node(370,2) = 191;
node(371,1) = 191; node(371,2) = 192;
node(372,1) = 192; node(372,2) = 193;
node(373,1) = 193; node(373,2) = 194;
node(374,1) = 194; node(374,2) = 195;
node(375,1) = 195; node(375,2) = 196;
node(376,1) = 196; node(376,2) = 197;
node(377,1) = 197; node(377,2) = 198;
node(378,1) = 198; node(378,2) = 199;
node(379,1) = 199; node(379,2) = 200;
node(380,1) = 200; node(380,2) = 201;
node(381,1) = 201; node(381,2) = 202;
node(382,1) = 202; node(382,2) = 203;
node(383,1) = 203; node(383,2) = 204;
node(384,1) = 204; node(384,2) = 205;
node(385,1) = 205; node(385,2) = 206;
node(386,1) = 206; node(386,2) = 207;
node(387,1) = 207; node(387,2) = 208;
node(388,1) = 208; node(388,2) = 209;
node(389,1) = 209; node(389,2) = 210;

node(390,1) = 190; node(390,2) = 211;
node(391,1) = 191; node(391,2) = 212;
node(392,1) = 192; node(392,2) = 213;
node(393,1) = 193; node(393,2) = 214;
node(394,1) = 194; node(394,2) = 215;
node(395,1) = 195; node(395,2) = 216;
node(396,1) = 196; node(396,2) = 217;
node(397,1) = 197; node(397,2) = 218;
node(398,1) = 198; node(398,2) = 219;
node(399,1) = 199; node(399,2) = 220;
node(400,1) = 200; node(400,2) = 221;
node(401,1) = 201; node(401,2) = 222;
node(402,1) = 202; node(402,2) = 223;
node(403,1) = 203; node(403,2) = 224;
node(404,1) = 204; node(404,2) = 225;
node(405,1) = 205; node(405,2) = 226;
node(406,1) = 206; node(406,2) = 227;
node(407,1) = 207; node(407,2) = 228;
node(408,1) = 208; node(408,2) = 229;
node(409,1) = 209; node(409,2) = 230;
node(410,1) = 210; node(410,2) = 231;

node(411,1) = 211; node(411,2) = 212;
node(412,1) = 212; node(412,2) = 213;
node(413,1) = 213; node(413,2) = 214;
node(414,1) = 214; node(414,2) = 215;
node(415,1) = 215; node(415,2) = 216;
node(416,1) = 216; node(416,2) = 217;
node(417,1) = 217; node(417,2) = 218;
node(418,1) = 218; node(418,2) = 219;
node(419,1) = 219; node(419,2) = 220;
node(420,1) = 220; node(420,2) = 221;
node(421,1) = 221; node(421,2) = 222;
node(422,1) = 222; node(422,2) = 223;
node(423,1) = 223; node(423,2) = 224;
node(424,1) = 224; node(424,2) = 225;
node(425,1) = 225; node(425,2) = 226;
node(426,1) = 226; node(426,2) = 227;
node(427,1) = 227; node(427,2) = 228;
node(428,1) = 228; node(428,2) = 229;
node(429,1) = 229; node(429,2) = 230;
node(430,1) = 230; node(430,2) = 231;

node(431,1) = 211; node(431,2) = 232;
node(432,1) = 212; node(432,2) = 233;
node(433,1) = 213; node(433,2) = 234;
node(434,1) = 214; node(434,2) = 235;
node(435,1) = 215; node(435,2) = 236;
node(436,1) = 216; node(436,2) = 237;
node(437,1) = 217; node(437,2) = 238;
node(438,1) = 218; node(438,2) = 239;
node(439,1) = 219; node(439,2) = 240;
node(440,1) = 220; node(440,2) = 241;
node(441,1) = 221; node(441,2) = 242;
node(442,1) = 222; node(442,2) = 243;
node(443,1) = 223; node(443,2) = 244;
node(444,1) = 224; node(444,2) = 245;
node(445,1) = 225; node(445,2) = 246;
node(446,1) = 226; node(446,2) = 247;
node(447,1) = 227; node(447,2) = 248;
node(448,1) = 228; node(448,2) = 249;
node(449,1) = 229; node(449,2) = 250;
node(450,1) = 230; node(450,2) = 251;
node(451,1) = 231; node(451,2) = 252;

node(452,1) = 232; node(452,2) = 233;
node(453,1) = 233; node(453,2) = 234;
node(454,1) = 234; node(454,2) = 235;
node(455,1) = 235; node(455,2) = 236;
node(456,1) = 236; node(456,2) = 237;
node(457,1) = 237; node(457,2) = 238;
node(458,1) = 238; node(458,2) = 239;
node(459,1) = 239; node(459,2) = 240;
node(460,1) = 240; node(460,2) = 241;
node(461,1) = 241; node(461,2) = 242;
node(462,1) = 242; node(462,2) = 243;
node(463,1) = 243; node(463,2) = 244;
node(464,1) = 244; node(464,2) = 245;
node(465,1) = 245; node(465,2) = 246;
node(466,1) = 246; node(466,2) = 247;
node(467,1) = 247; node(467,2) = 248;
node(468,1) = 248; node(468,2) = 249;
node(469,1) = 249; node(469,2) = 250;
node(470,1) = 250; node(470,2) = 251;
node(471,1) = 251; node(471,2) = 252;

node(472,1) = 232; node(472,2) = 253;
node(473,1) = 233; node(473,2) = 254;
node(474,1) = 234; node(474,2) = 255;
node(475,1) = 235; node(475,2) = 256;
node(476,1) = 236; node(476,2) = 257;
node(477,1) = 237; node(477,2) = 258; 
node(478,1) = 238; node(478,2) = 259;
node(479,1) = 239; node(479,2) = 260;
node(480,1) = 240; node(480,2) = 261;
node(481,1) = 241; node(481,2) = 262;
node(482,1) = 242; node(482,2) = 263;
node(483,1) = 243; node(483,2) = 264;
node(484,1) = 244; node(484,2) = 265;
node(485,1) = 245; node(485,2) = 266;
node(486,1) = 246; node(486,2) = 267;
node(487,1) = 247; node(487,2) = 268;
node(488,1) = 248; node(488,2) = 269;
node(489,1) = 249; node(489,2) = 270;
node(490,1) = 250; node(490,2) = 271;
node(491,1) = 251; node(491,2) = 272;
node(492,1) = 252; node(492,2) = 273;

node(493,1) = 253; node(493,2) = 254;
node(494,1) = 254; node(494,2) = 255;
node(495,1) = 255; node(495,2) = 256;
node(496,1) = 256; node(496,2) = 257;
node(497,1) = 257; node(497,2) = 258;
node(498,1) = 258; node(498,2) = 259;
node(499,1) = 259; node(499,2) = 260;
node(500,1) = 260; node(500,2) = 261;
node(501,1) = 261; node(501,2) = 262;
node(502,1) = 262; node(502,2) = 263;
node(503,1) = 263; node(503,2) = 264;
node(504,1) = 264; node(504,2) = 265;
node(505,1) = 265; node(505,2) = 266;
node(506,1) = 266; node(506,2) = 267;
node(507,1) = 267; node(507,2) = 268;
node(508,1) = 268; node(508,2) = 269;
node(509,1) = 269; node(509,2) = 270;
node(510,1) = 270; node(510,2) = 271;
node(511,1) = 271; node(511,2) = 272;
node(512,1) = 272; node(512,2) = 273;

node(513,1) = 253; node(513,2) = 274;
node(514,1) = 254; node(514,2) = 275;
node(515,1) = 255; node(515,2) = 276;
node(516,1) = 256; node(516,2) = 277;
node(517,1) = 257; node(517,2) = 278;
node(518,1) = 258; node(518,2) = 279;
node(519,1) = 259; node(519,2) = 280;
node(520,1) = 260; node(520,2) = 281;
node(521,1) = 261; node(521,2) = 282;
node(522,1) = 262; node(522,2) = 283;
node(523,1) = 263; node(523,2) = 284;
node(524,1) = 264; node(524,2) = 285;
node(525,1) = 265; node(525,2) = 286;
node(526,1) = 266; node(526,2) = 287;
node(527,1) = 267; node(527,2) = 288;
node(528,1) = 268; node(528,2) = 289;
node(529,1) = 269; node(529,2) = 290;
node(530,1) = 270; node(530,2) = 291;
node(531,1) = 271; node(531,2) = 292;
node(532,1) = 272; node(532,2) = 293;
node(533,1) = 273; node(533,2) = 294;

node(534,1) = 274; node(534,2) = 275;
node(535,1) = 275; node(535,2) = 276;
node(536,1) = 276; node(536,2) = 277;
node(537,1) = 277; node(537,2) = 278;
node(538,1) = 278; node(538,2) = 279;
node(539,1) = 279; node(539,2) = 280;
node(540,1) = 280; node(540,2) = 281;
node(541,1) = 281; node(541,2) = 282;
node(542,1) = 282; node(542,2) = 283;
node(543,1) = 283; node(543,2) = 284;
node(544,1) = 284; node(544,2) = 285;
node(545,1) = 285; node(545,2) = 286;
node(546,1) = 286; node(546,2) = 287;
node(547,1) = 287; node(547,2) = 288;
node(548,1) = 288; node(548,2) = 289;
node(549,1) = 289; node(549,2) = 290;
node(550,1) = 290; node(550,2) = 291;
node(551,1) = 291; node(551,2) = 292;
node(552,1) = 292; node(552,2) = 293;
node(553,1) = 293; node(553,2) = 294;

node(554,1) = 274; node(554,2) = 295;
node(555,1) = 275; node(555,2) = 296;
node(556,1) = 276; node(556,2) = 297;
node(557,1) = 277; node(557,2) = 298;
node(558,1) = 278; node(558,2) = 299;
node(559,1) = 279; node(559,2) = 300;
node(560,1) = 280; node(560,2) = 301;
node(561,1) = 281; node(561,2) = 302;
node(562,1) = 282; node(562,2) = 303;
node(563,1) = 283; node(563,2) = 304;
node(564,1) = 284; node(564,2) = 305;
node(565,1) = 285; node(565,2) = 306;
node(566,1) = 286; node(566,2) = 307;
node(567,1) = 287; node(567,2) = 308;
node(568,1) = 288; node(568,2) = 309;
node(569,1) = 289; node(569,2) = 310;
node(570,1) = 290; node(570,2) = 311;
node(571,1) = 291; node(571,2) = 312;
node(572,1) = 292; node(572,2) = 313;
node(573,1) = 293; node(573,2) = 314;
node(574,1) = 294; node(574,2) = 315;

node(575,1) = 295; node(575,2) = 296;
node(576,1) = 296; node(576,2) = 297;
node(577,1) = 297; node(577,2) = 298;
node(578,1) = 298; node(578,2) = 299;
node(579,1) = 299; node(579,2) = 300;
node(580,1) = 300; node(580,2) = 301;
node(581,1) = 301; node(581,2) = 302;
node(582,1) = 302; node(582,2) = 303;
node(583,1) = 303; node(583,2) = 304;
node(584,1) = 304; node(584,2) = 305;
node(585,1) = 305; node(585,2) = 306;
node(586,1) = 306; node(586,2) = 307;
node(587,1) = 307; node(587,2) = 308;
node(588,1) = 308; node(588,2) = 309;
node(589,1) = 309; node(589,2) = 310;
node(590,1) = 310; node(590,2) = 311;
node(591,1) = 311; node(591,2) = 312;
node(592,1) = 312; node(592,2) = 313;
node(593,1) = 313; node(593,2) = 314;
node(594,1) = 314; node(594,2) = 315;

node(595,1) = 295; node(595,2) = 316;
node(596,1) = 296; node(596,2) = 317;
node(597,1) = 297; node(597,2) = 318;
node(598,1) = 298; node(598,2) = 319;
node(599,1) = 299; node(599,2) = 320;
node(600,1) = 300; node(600,2) = 321;
node(601,1) = 301; node(601,2) = 322;
node(602,1) = 302; node(602,2) = 323;
node(603,1) = 303; node(603,2) = 324;
node(604,1) = 304; node(604,2) = 325;
node(605,1) = 305; node(605,2) = 326;
node(606,1) = 306; node(606,2) = 327;
node(607,1) = 307; node(607,2) = 328;
node(608,1) = 308; node(608,2) = 329;
node(609,1) = 309; node(609,2) = 330;
node(610,1) = 310; node(610,2) = 331;
node(611,1) = 311; node(611,2) = 332;
node(612,1) = 312; node(612,2) = 333;
node(613,1) = 313; node(613,2) = 334;
node(614,1) = 314; node(614,2) = 335;
node(615,1) = 315; node(615,2) = 336;

node(616,1) = 316; node(616,2) = 317;
node(617,1) = 317; node(617,2) = 318;
node(618,1) = 318; node(618,2) = 319;
node(619,1) = 319; node(619,2) = 320;
node(620,1) = 320; node(620,2) = 321;
node(621,1) = 321; node(621,2) = 322;
node(622,1) = 322; node(622,2) = 323;
node(623,1) = 323; node(623,2) = 324;
node(624,1) = 324; node(624,2) = 325;
node(625,1) = 325; node(625,2) = 326;
node(626,1) = 326; node(626,2) = 327;
node(627,1) = 327; node(627,2) = 328;
node(628,1) = 328; node(628,2) = 329;
node(629,1) = 329; node(629,2) = 330;
node(630,1) = 330; node(630,2) = 331;
node(631,1) = 331; node(631,2) = 332;
node(632,1) = 332; node(632,2) = 333;
node(633,1) = 333; node(633,2) = 334;
node(634,1) = 334; node(634,2) = 335;
node(635,1) = 335; node(635,2) = 336;

node(636,1) = 316; node(636,2) = 337;
node(637,1) = 317; node(637,2) = 338;
node(638,1) = 318; node(638,2) = 339;
node(639,1) = 319; node(639,2) = 340;
node(640,1) = 320; node(640,2) = 341;
node(641,1) = 321; node(641,2) = 342;
node(642,1) = 322; node(642,2) = 343;
node(643,1) = 323; node(643,2) = 344;
node(644,1) = 324; node(644,2) = 345;
node(645,1) = 325; node(645,2) = 346;
node(646,1) = 326; node(646,2) = 347;
node(647,1) = 327; node(647,2) = 348;
node(648,1) = 328; node(648,2) = 349;
node(649,1) = 329; node(649,2) = 350;
node(650,1) = 330; node(650,2) = 351;
node(651,1) = 331; node(651,2) = 352;
node(652,1) = 332; node(652,2) = 353;
node(653,1) = 333; node(653,2) = 354;
node(654,1) = 334; node(654,2) = 355;
node(655,1) = 335; node(655,2) = 356;
node(656,1) = 336; node(656,2) = 357;

node(657,1) = 337; node(657,2) = 338;
node(658,1) = 338; node(658,2) = 339;
node(659,1) = 339; node(659,2) = 340;
node(660,1) = 340; node(660,2) = 341;
node(661,1) = 341; node(661,2) = 342;
node(662,1) = 342; node(662,2) = 343;
node(663,1) = 343; node(663,2) = 344;
node(664,1) = 344; node(664,2) = 345;
node(665,1) = 345; node(665,2) = 346;
node(666,1) = 346; node(666,2) = 347;
node(667,1) = 347; node(667,2) = 348;
node(668,1) = 348; node(668,2) = 349;
node(669,1) = 349; node(669,2) = 350;
node(670,1) = 350; node(670,2) = 351;
node(671,1) = 351; node(671,2) = 352;
node(672,1) = 352; node(672,2) = 353;
node(673,1) = 353; node(673,2) = 354;
node(674,1) = 354; node(674,2) = 355;
node(675,1) = 355; node(675,2) = 356;
node(676,1) = 356; node(676,2) = 357;

node(677,1) = 337; node(677,2) = 358;
node(678,1) = 338; node(678,2) = 359;
node(679,1) = 339; node(679,2) = 360;
node(680,1) = 340; node(680,2) = 361;
node(681,1) = 341; node(681,2) = 362;
node(682,1) = 342; node(682,2) = 363;
node(683,1) = 343; node(683,2) = 364;
node(684,1) = 344; node(684,2) = 365;
node(685,1) = 345; node(685,2) = 366;
node(686,1) = 346; node(686,2) = 367;
node(687,1) = 347; node(687,2) = 368;
node(688,1) = 348; node(688,2) = 369;
node(689,1) = 349; node(689,2) = 370;
node(690,1) = 350; node(690,2) = 371;
node(691,1) = 351; node(691,2) = 372;
node(692,1) = 352; node(692,2) = 373;
node(693,1) = 353; node(693,2) = 374;
node(694,1) = 354; node(694,2) = 375;
node(695,1) = 355; node(695,2) = 376;
node(696,1) = 356; node(696,2) = 377;
node(697,1) = 357; node(697,2) = 378;

for i = 358:377
    node(i + 340, 1) = i;
    node(i + 340, 2) = i + 1;
end

for i = 718:738
    node(i, 1) = i - 360;
    node(i, 2) = i - 339;
end

for i = 739:758
    node(i, 1) = i - 360;
    node(i, 2) = i - 359;
end

for i = 759:779
	node(i, 1) = i - 380;
    node(i, 2) = i - 359;
end

for i = 780:799
    node(i, 1) = i - 380;
    node(i, 2) = i - 379;
end

for i = 800:820
	node(i, 1) = i - 400;
    node(i, 2) = i - 379;
end

for i = 821:840
    node(i, 1) = i - 400;
    node(i, 2) = i - 399;
end  

nbc=12;  % So luong chuyen vi nut bi rang buoc
% Nut 1. Chuyen vi nut bi rang buoc la q1, q3, q6
nq1 = 1;
nq2 = 2;
nq3 = 3;
i = 1;
bcdof(1)=(i-1)*6 + nq1; bcval(1)=0;
bcdof(2)=(i-1)*6 + nq2; bcval(2)=0;
bcdof(3)=(i-1)*6 + nq3; bcval(3)=0;
% Nut 3
i = 21;
bcdof(4)=(i-1)*6 + nq1; bcval(4)=0;
bcdof(5)=(i-1)*6 + nq2; bcval(5)=0;
bcdof(6)=(i-1)*6 + nq3; bcval(6)=0;
% Nut 7
i = 421;
bcdof(7)=(i-1)*6 + nq1; bcval(7)=0;
bcdof(8)=(i-1)*6 + nq2; bcval(8)=0;
bcdof(9)=(i-1)*6 + nq3; bcval(9)=0;
% Nut 9
i = 441;
bcdof(10)=(i-1)*6 + nq1; bcval(10)=0;
bcdof(11)=(i-1)*6 + nq2; bcval(11)=0;
bcdof(12)=(i-1)*6 + nq3; bcval(12)=0;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Bac tu do, tong so nut, tong so bac tu do, tong so doan thanh (phan tu)
% q1 - Chuyen vi doc truc tai nut 1 (u_1). *
% q2 - Chuyen vi doc truc tai nut 2 (u_2).
% q3 - Chuyen vi ngang theo phuong y tai nut 1 (v_1). *
% q4 - Chuyen vi ngang theo phuong z tai nut 2 (v_2).
% q5 - Chuyen vi ngang theo phuong z tai nut 1 (w_1).
% q6 - Chuyen vi ngang theo phuong y tai nut 2 (w_2). *
% q7 - Goc xoan tai nut 1 (teta_1x).
% q8 - Goc xoan tai nut 2 (teta_2x).
% q9 - Goc xoay quan truc y tai dau 1 (teta_1y).
% q10 - Goc xoay quan truc y tai dau 2 (teta_2y).
% q11 - Goc xoay quan truc z tai dau 1 (teta_1z).
% q12 - Goc xoay quan truc z tai dau 2 (teta_2z).

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Thanh phan noi luc tai 2 nut cua (phan tu) doan thanh.
% N1 - luc doc tai nut 1 (N_1).
% N2 - luc doc tai nut 2 (N_2).
% Q1y - luc cat tren phuong y tai nut 1.
% Q2y - luc cat tren phuong y tai nut 2.
% Q1z - luc cat tren phuong z tai nut 1.
% Q2z - luc cat tren phuong z tai nut 2.
% M1x - mo men xoan tai nut 1.
% M2x - mo men xoan tai nut 2.
% M1y - mo men quanh truc y tai nut 1.
% M2y - mo men quanh truc y tai nut 2.
% M1z - mo men quanh truc z tai nut 1.
% M2z - mo men quanh truc z tai nut 2.

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Xay dung ma tran do cung cua thanh hai dau han trong khong gian 3 chieu
% Ma tran do cung dia phuong (cuc bo) co ma tran 12 x 12
S1 = E*Area/L;
S2 = 12*E*Iz/L^3;
S3 = 6*E*Iz/L^2;
S4 = 12*E*Iy/L^3;
S5 = 6*E*Iy/L^2;
S6 = G*Ix/L;
S7 = 4*E*Iy/L;
S8 = 2*E*Iy/L;
S9 = 4*E*Iz/L;
S10 = 2*E*Iz/L;

k_local = [S1, 0, 0, 0, 0, 0, -S1, 0, 0, 0, 0, 0; ... 
           0, S2, 0, 0, 0, S3, 0, -S2, 0, 0, 0, S3; ... 
           0, 0, S4, 0, -S5, 0, 0, 0, -S4, 0, -S5, 0; ... 
           0, 0, 0, S6, 0, 0, 0, 0, 0, -S6, 0, 0; ... 
           0, 0, -S5, 0, S7, 0, 0, 0, S5, 0, S8, 0; ... 
           0, S3, 0, 0, 0, S9, 0, -S3, 0, 0, 0, S10; ... 
           -S1, 0, 0, 0, 0, 0, S1, 0, 0, 0, 0, 0; ... 
           0, -S2, 0, 0, 0, -S3, 0, S2, 0, 0, 0, -S3; ... 
           0, 0, -S4, 0, S5, 0, 0, 0, S4, 0, S5, 0; ... 
           0, 0, 0, -S6, 0, 0, 0, 0, 0, S6, 0, 0; ... 
           0, 0, -S5, 0, S8, 0, 0, 0, S5, 0, S7, 0; ... 
           0, S3, 0, 0, 0, S10, 0, -S3, 0, 0, 0, S9]; 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Ghep noi ma tran do cung dia phuong thanh ma tran do cung tong the

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Gan gia tri 0
ff=zeros(sdof,1);
kk=zeros(sdof,sdof);
rf=zeros(sdof,1);
index=zeros(nnel*ndof,1);
noiluc=zeros(12,nel);
cvnutpt=zeros(12,1);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Khai bao vec to tai quy nut tai cac diem giao dam
% co the bang phuong phap ngau nhien - xac suat thong ke
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Tai trong quy nut dam
q = 500 * 0.6 * 0.6; % Tai trong theo phuong dung tai nut so 5
for i = 1:441
    ff((i-1)*ndof+nq2) = -0.18;
%     ff((i-1)*ndof+nq2) = -rand(1);
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%Xac dinh ma tran do cung cua he
for iel=1:nel
    index=mahoacvinut(iel,node,nnel,ndof);
    
    nd(1)=node(iel,1);
    nd(2)=node(iel,2);
    
    x1=x(nd(1));z1=z(nd(1));y1=y(nd(1));
    x2=x(nd(2));z2=z(nd(2));y2=y(nd(2));
    
    leng=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    if (x2-x1)==0
        beta=pi/2;
    else
        beta=atan((y2-y1)/(x2-x1));
    end

%     L = ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)^0.5

    % Chuyen ma tran tu he toa do dia phong ve he toa do tong the

Cx = (x2 - x1)/L;
Cy = (y2 - y1)/L;
Cz = (z2 - z1)/L;

cos_beta = cos(0);
sin_beta = sin(0);

cos_gamma = (Cx^2 + Cz^2);
sin_gamma = Cy;

alpha = 0;

T_alpha = [1,0,0;...
    0,cos(alpha),sin(alpha);...
    0,-sin(alpha),cos(alpha)];

T_gamma = [cos_gamma, sin_gamma, 0;...
    -sin_gamma, cos_gamma, 0;...
    0,0,1];

T_beta = [1, 0, 0; ...
    0,1,0;...
    0, 0, 1];

zeros_T = zeros(3,3);

T = T_alpha * T_gamma * T_beta;
% 
% % r = horzcat(T,  zeros_T,    zeros_T,    zeros_T);
r = [T,   zeros_T,    zeros_T,    zeros_T; ...
    zeros_T,    T,  zeros_T,    zeros_T; ...
    zeros_T,    zeros_T,    T,  zeros_T; ...
    zeros_T,    zeros_T,    zeros_T,    T];

%     r=mtctoadoUBoot(alpha, beta, gamma);
    [k_local] = mtdcdam(E,Area,Iz,Iy,Ix,L);
    k = k_local;
%     [k]=mtdcdamUBoot(el,area,xi,b,ko,leng);
    k=r'*k*r;
    [kk]=mtdche(kk,k,index);
%     [kk]=mtdche(kk,k,index);
end
% Khu dieu kien bien
for i = 1:sdof
    ff(i)=ff(i)+rf(i);
end
[kk,ff]=khudkbien(kk,ff,bcdof,bcval);
% [kk,ff]=khudkbien(kk,ff,bcdof,bcval)

disp=inv(kk)*ff;
for iel=1:nel
    nd(1)=node(iel,1);
    nd(2)=node(iel,2);
    x1=x(nd(1));z1=z(nd(1));
    x2=x(nd(2));z2=z(nd(2));
    [k_local] = mtdcdam(E,Area,Iz,Iy,Ix,L);
    k = k_local;
%     k=mtdcdamnen(el,area,xi,b,ko,leng);
    % cvi nut PT
    cvinutpt(1,1)=disp(6*nd(1)-5,1);
    cvinutpt(2,1)=disp(6*nd(1)-4,1);
    cvinutpt(3,1)=disp(6*nd(1)-3,1);
    cvinutpt(4,1)=disp(6*nd(1)-2,1);
    cvinutpt(5,1)=disp(6*nd(1)-1,1);
    cvinutpt(6,1)=disp(6*nd(1)-0,1);
    cvinutpt(7,1)=disp(6*nd(2)-5,1);
    cvinutpt(8,1)=disp(6*nd(2)-4,1);
    cvinutpt(9,1)=disp(6*nd(2)-3,1);
    cvinutpt(10,1)=disp(6*nd(2)-2,1);
    cvinutpt(11,1)=disp(6*nd(2)-1,1);
    cvinutpt(12,1)=disp(6*nd(2)-0,1);
    cvinutpt=r*cvinutpt; % chuyen ve Hetoado diaphuong
    % Noi luc PT khong tinh tai tac dung trong Phantu
    % Noi luc mo men 8 va 12
    noiluc(:,iel)=k*cvinutpt;
end

A = kk; b = ff; A = sparse(A);
[a_val, a_row_ptr, a_col_idx] = sparse2csr(A,1);
% Chon tien dieu kien phai K
p = 13; a = kk;
a = ILU_P(a, p);
% K = l * l';
K = a;
[k_val, k_row_ptr, k_col_idx] = sparse2csr(K,1);
% Khoi tao gia tri x_0 ban dau
n = size(A,1);
x_0 = zeros(n,1);
% r_0 = b - A* x_0;
r_0 = b - dotproduct(a_val, a_row_ptr, a_col_idx, x_0);
r_bar = r_0;
i = 1;
tol = 1e-10;
rho_0 = r_bar' * r_0;
p_1 = r_0;
% Giai he phuong trinh K*p_hat = p_1 theo phuong phap GMRES (tam thoi dung)
% res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
% m_diagonal = eye(n,1); res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
m_diagonal = diag(A); res_tol = 1e-10; max_iter = 1000;  use_HH = 0; restart = 10;
% [p_hat, ~, ~, ~] = GMRES(K, p_1, restart, res_tol, max_iter, M, use_HH);
[p_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, p_1, restart, res_tol, max_iter, m_diagonal, use_HH);

% v_1 = A * p_hat;
v_1 = dotproduct(a_val, a_row_ptr, a_col_idx, p_hat);
alpha_1 = rho_0 / (r_bar' * v_1);
s = r_0 - alpha_1 * v_1;
% [s_hat, ~, ~, ~] = GMRES(K, s, restart, res_tol, max_iter, M, use_HH);
[s_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, s, restart, res_tol, max_iter, m_diagonal, use_HH);

% t = A * s_hat;
t = dotproduct(a_val, a_row_ptr, a_col_idx, s_hat);

omega_1 = t' *s / (t' * t);
x_1 =  x_0 + alpha_1 * p_hat + omega_1 * s_hat;
r_1 = s - omega_1 * t;
r_bar = r_1;
i = i + 1;
res = x_1 - x_0;
norm_res = norm(res,2);
while ((i < max_iter) && (norm_res > tol))
    rho_1 = r_bar' * r_0;
    beta_1 = (rho_1 / rho_0) * (alpha_1 / omega_1);
    p_2 = r_1 + beta_1 * (p_1 - omega_1 * v_1);
%     [p_hat, ~, ~, ~] = GMRES(K, p_2, restart, res_tol, max_iter, M, use_HH);
    [p_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, p_2, restart, res_tol, max_iter, m_diagonal, use_HH);
    
%     v_2 = A * p_hat;
    v_2 = dotproduct(a_val, a_row_ptr, a_col_idx, p_hat);
    
    alpha_2 = rho_1 / (r_bar' * v_2);
    s = r_1 - alpha_2 * v_2;
%     [s_hat, ~, ~, ~] = GMRES(K, s, restart, res_tol, max_iter, M, use_HH);
    [s_hat, ~, ~, ~] = GMRES_Diagonal(k_val, k_row_ptr, k_col_idx, s, restart, res_tol, max_iter, m_diagonal, use_HH);

%     t = A * s_hat;
    t = dotproduct(a_val, a_row_ptr, a_col_idx, s_hat);
    
    omega_2 = t' * s / (t' * t);
    x_2 = x_1 + alpha_2 * p_hat + omega_2 * s_hat;
    r_2 = s - omega_2 * t;
    r_1 = r_2;
    x_1 = x_2;
    p_1 = p_2;
    rho_0 = rho_1;
    alpha_1 = alpha_2;
    omega_1 = omega_2;
    v_1 = v_2;
    r_bar = r_1;
    i = i + 1;
    res = x_2 - x_1;
    norm_res = norm(res,2);
end
rr = A * x_2 - b;
norm(rr,2)