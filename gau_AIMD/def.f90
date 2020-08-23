! This file is a part of wQC(https://github.com/wubaihua/wQC)

! Author:
! > Baihua Wu
! > wubaihua@pku.edu.cn

! Last update: 2020-5-3

! module def: some definition in wQC, including the type of atom
! and shell, degeneracy of shell in different type of GTO, etc. 

! Some code in this file refers to the source code of Multiwfn 3.6
! (http://sobereva.com/multiwfn/). Multiwfn is a quantum chemistry
! wavefunction analysis program developed by Dr. Tian Lu. Author 
! thanks for his contribution.



module def
    type atomtype !�Զ���ļ�¼ԭ����Ϣ������
    character*2 name !ԭ����
    integer index !ԭ����š���ʹ��������ʱԭ�Ӻ˵�ɽ�С��ԭ�����
    real*8 x,y,z,charge !ԭ������(Bohr)��ԭ�Ӻ˵��
    end type

    type primtype !�Զ���ļ�¼GTF��Ϣ��ר������
    integer center,functype !GTF����ԭ���Լ�GTF������
    real*8,allocatable:: e(:)!ָ��

    real*8,allocatable::  C(:,:)
    end type

character*2 :: name2ind(0:109)=(/ "Bq","H ","He", &   !ԭ������ԭ����ŵ�ת������0������ԭ��
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt" /) !104~109


character*2 :: shellindex(1:6)=(/'S ','P ','SP','D','F','G'/)

!����������GTF������GTF��x,y,z�ϵ�ָ����ת����������XY��GTF��x,y,z��ָ���ֱ���1,1,0�����������еĵ�8�ж�Ӧ������������GTF��functype��8��
integer :: type2ix(35)=(/ 0,1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4 /)
integer :: type2iy(35)=(/ 0,0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
integer :: type2iz(35)=(/ 0,0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
integer :: nmo=0,nprims=0,ncenter=0 !�������GTF����ԭ����
integer wfntype !0/1/2������������R/U/ROHF��������3/4�ֱ�����տǲ�Ϳ��ǲ��HF������
real*8 :: nelec=0,naelec=0,nbelec=0 !�ܵ�������alpha��beta������
type(atomtype),allocatable :: a(:) !��¼ԭ����Ϣ������
type(primtype),allocatable :: b(:) !��¼GTF��Ϣ������
real*8,allocatable :: MOocc(:),MOene(:) !���ռ�����͹������
integer,allocatable :: MOtype(:) !��¼������͡�0/1/2�ֱ����������/alpha/beta�͹��
!real*8,allocatable :: CO(:,:) !���չ��ϵ������CO(i,j)������j��GTF�ڵ�i�Ź���е�չ��ϵ����ϵ�����Ѿ�������ϵ������һ��ϵ��ȫ��������ȥ�ˡ�
!index of shell : 0--S,  1--P,  2--D,  3--F     
integer :: shl_degen_sph(0:3)=(/1,3,5,7/) !Degeneracy of shell in Spherical harmonic Gaussian function(5D,7F) 
integer :: shl_degen_cart(0:3)=(/1,3,6,10/) !Degeneracy of shell in Cartesian Gaussian function(6D,10F) 

contains


 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end module


