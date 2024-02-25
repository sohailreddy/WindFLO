!
!
!	Copyright 2019, Sohail R. Reddy
!	email: sredd001@fiu.edu
!	www.sohailreddy.com
!
!

module WindFLO

	use iso_c_binding
	use octree_mod
	use terrain_mod
	
	implicit none

	real(8), parameter :: undefined = 9876543210.0d0
	integer, parameter :: undefined_i = 987654321
	real(8), parameter :: pi = 3.14159265359d0
	real(8), parameter :: vonKarmanConstant = 0.40d0
	real(8) :: viscosity = 1.81d-5



type windRoseProps
	real(8) :: modelVelocity(3)
	real(8) :: probability
end type windRoseProps


type turbineProps

	character(len = 100) :: name ! name of turbine
	integer :: turbineNum		! the ith turbine
	real(8) :: position(3)		! x, y, z position of turbine
	real(8) :: orientation(3)	! polar orientation of turbine	
	real(8) :: height			! rotor height
	real(8) :: radius			! radius of rotor
	real(8) :: diameter 		! diameter of rotor
	real(8) :: ratedPower		! rated power of the turbine

!	misc
	logical :: fictitious 		! is it an imaginary turbine
	logical :: yaw						! is the turbine able to yaw
	
	real(8), allocatable :: CpCurve(:,:)			! Cp curve data for spline fitting
	real(8), allocatable :: CtCurve(:,:)			! Ct curve data for spline fitting

end type turbineProps




type WindFLO_params

	character(len = 1000) :: fileName;
	integer nTurbines;
	logical batch;

!	Atmospheric parameters
	real(8) rho;
	real(8) turbulenceIntensity;

!	Ambient wind parameters	
	character(len = 100) :: windModel;	! wind model
	real(8) :: modelVelocity(3);		! can be for constant, friction velocity (for log) or ref velocity (for power)
	real(8) :: surfaceRoughness;		! for log profile
	real(8) :: referenceHeight;			! for power profile

!	Wake modeling parameters
	character(len = 100) wakeModel;			! wake model
	character(len = 100) wakeMergeModel;	! wake merge model
	real(8) :: wakeExpansionCoeff(2);		! wake expansion coefficient ... 2 for XA model
	
!	Numerical parameters
	integer gaussOrder;
	integer monteCarloPts;
	
!	Turbine parameters
	type(TurbineProps), allocatable :: turbines(:)
	
!	Cost parameters
	real(8) :: coe	! Cost of energy

!	Terrain parameters
	character(len = 1000) terrainModel	! Terrain interpolation model	
	
!	Rotation Matrix
	real(8) :: rotMat(3,3)
	real(8) :: rotMat_inv(3,3)

!	Wind Rose Data
	type(windRoseProps), allocatable :: windRose(:)
	
	logical :: alreadyRead = .false.
	
	
end type WindFLO_params





!	Analysis params
	integer totTurbines						! total number of turbines
	character(len=1000) :: offBodyPoints	! off body points at which to compute velocities if needed


	type(turbineProps), allocatable :: turbines(:)
	type(WindFLO_params) :: WindFLO_input
	real(8), allocatable :: terrain(:,:)

	character(len = 1000) :: terrainFileCurrent = ''
	character(len = 1000) :: windRoseFileCurrent = ''




contains


!------------------------------------------------------------
!
!	Preprocessing and Read Input
!
!------------------------------------------------------------





subroutine ReadInput(inFileName) bind(C,name='ReadInput_')

	implicit none
	character(c_char) :: inFileName(1000)
	character(len = 1000) :: fileName = ''
	integer i

	integer nTurbines;
	logical batch;

!	Atmospheric parameters
	real(8) rho;
	real(8) turbulenceIntensity;

!	Ambient wind parameters	
	character(len = 100) :: windModel;	! wind model
	real(8) :: modelVelocity(3);		! can be for constant, friction velocity (for log) or ref velocity (for power)
	real(8) :: surfaceRoughness;		! for log profile
	real(8) :: referenceHeight;			! for power profile

!	Wake modeling parameters
	character(len = 100) wakeModel;			! wake model
	character(len = 100) wakeMergeModel;	! wake merge model
	real(8) :: wakeExpansionCoeff(2);		! wake expansion coefficient ... 2 for XA model
	
!	Numerical parameters
	integer gaussOrder;
	integer monteCarloPts;

!	Cost parameters
	real(8) :: coe	

!	Terrain parameters
!	integer :: nNearest		NOT USED ! Already defined in terrain_mod
!	integer :: rbfKernel	! Already defined in terrain_mod
!	real(8) :: shapeFactor	! Already defined in terrain_mod
	character(len = 1000) terrainModel	! Terrain interpolation model
	
	character(len = 1000) :: turbineFiles(1000) = ''
	character(len = 1000) :: terrainFile = ''
	character(len = 1000) :: windRoseFile = ''

	!------------------------------------------------
	namelist / 	WindFLO_data / &
				nTurbines, batch,								&
				rho, turbulenceIntensity,						&
				windModel, modelVelocity, 						&
				surfaceRoughness, referenceHeight,				&
				wakeModel, wakeMergeModel, wakeExpansionCoeff,	&
				gaussOrder, monteCarloPts, 						&
				turbineFiles,									&
				coe,											&
				nNearest, rbfKernel, shapeFactor, powerIDW,		&
				terrainModel, terrainFile, windRoseFile,					&
!	Misc Parameters
				octreeDepth, octreeMaxPts
	!------------------------------------------------				

!	if(WindFLO_input%alreadyRead) then
!		return
!	end if


	fileName = getFString(inFileName)
	open(unit=100, file = trim(fileName), status = 'unknown')				
	read(100,nml=WindFLO_data)
	close(100)

	nTurbines = count(turbineFiles .ne. '')

	WindFLO_input%nTurbines = nTurbines
	WindFLO_input%batch = batch
	WindFLO_input%rho = rho
	WindFLO_input%turbulenceIntensity = turbulenceIntensity
	WindFLO_input%windModel = windModel
	WindFLO_input%modelVelocity = modelVelocity
	WindFLO_input%surfaceRoughness = surfaceRoughness
	WindFLO_input%referenceHeight = referenceHeight
	WindFLO_input%wakeModel = wakeModel
	WindFLO_input%wakeMergeModel = wakeMergeModel
	WindFLO_input%wakeExpansionCoeff = 	wakeExpansionCoeff
	WindFLO_input%gaussOrder = gaussOrder
	WindFLO_input%monteCarloPts = monteCarloPts
	WindFLO_input%coe = coe	
		

	if(.not. allocated(WindFLO_input%turbines)) then
		allocate(WindFLO_input%turbines(nTurbines))
	end if
		
	do i = 1, nTurbines
		call readTurbineData(trim(turbineFiles(i)), WindFLO_input%turbines(i))
	end do

	WindFLO_input%terrainModel = terrainModel
	if(trim(terrainFile) /= '')  call readTerrain(terrainFile)

	if(trim(windRoseFile) /= '') call readWindRose(windRoseFile)
				
	WindFLO_input%alreadyRead = .true.
	return
end subroutine ReadInput



subroutine readTurbineData(inputFile, inTurbine)
	implicit none
	
	character(*) :: inputFile
	type(turbineProps) :: inTurbine
	integer nCount, i, j, ierr
						
						
	character(len = 30) :: name 		! name of turbine
	integer :: turbineNum = undefined_i	! the ith turbine
	real(8) :: position(3)				! x, y, z position of turbine
	real(8) :: orientation(3)			! polar orientation of turbine	
	real(8) :: height 					! rotor height
	real(8) :: radius 					! radius of rotor
	real(8) :: diameter 				! diameter of rotor
	real(8) :: ratedPower				! rated power of the turbine
	real(8) :: CpCurve(2,5000)			! Cp curve data for spline fitting
	real(8) :: CtCurve(2,5000)			! Ct curve data for spline fitting
	logical :: fictitious				! is the turbine imaginary
	logical :: yaw						! is the turbine able to yaw
	
	!------------------------------------------------
	namelist / turbine_data /&
			turbineNum,	name, fictitious,		&
			position, height, diameter, radius,	&
			orientation, yaw, ratedPower,		&
			CpCurve, CtCurve	
	!------------------------------------------------
	
		
	if(trim(inputFile) == '') return
	
	turbineNum = undefined_i	! the ith turbine
	name = ''
	position = undefined	! x, y, z position of turbine
	height = undefined		! rotor height
	radius = undefined		! radius of rotor
	diameter = undefined	! diameter of rotor
	ratedPower = -undefined
	CpCurve = undefined
	CtCurve = undefined
	fictitious = .false.
	yaw = .true. 
	orientation = 0.0d0; orientation(1) = -1.0d0
		
	call initializeTurbineProps(inTurbine)
		
	open(unit=101, file = trim(inputFile), status = 'unknown')				
	read(101,nml=turbine_data)
	close(101)
	
	inTurbine%turbineNum = turbineNum
	inTurbine%name = trim(name)
	inTurbine%position = position
	inTurbine%orientation = orientation/norm2(orientation)
	inTurbine%height = height
	inTurbine%radius = radius
	inTurbine%diameter = diameter
	inTurbine%ratedPower = ratedPower
	inTurbine%fictitious = fictitious
	inTurbine%yaw = yaw

	if(inTurbine%radius == undefined) inTurbine%radius = inTurbine%diameter/2.0d0
	if(inTurbine%diameter == undefined) inTurbine%diameter = inTurbine%radius*2.0d0		
	
	call getCpOrCtCurve(CpCurve,CtCurve)
	
	nCount = count(CpCurve /= undefined) / 2
	if(nCount /= 0) then
		allocate(inTurbine%CpCurve(2, nCount))
		inTurbine%CpCurve = CpCurve(: , 1 : nCount)
	end if
	
	nCount = count(CtCurve /= undefined) / 2		
	if(nCount /= 0) then
		allocate(inTurbine%CtCurve(2, nCount))		
		inTurbine%CtCurve = CtCurve(: , 1 : nCount)
	end if	
				


!	open(unit = 1001, file = 'cpct.dat', status = 'unknown')
!	do i = 1, ncount
!		write(1001,*) inTurbine%CtCurve(1, i), inTurbine%CpCurve(2, i),inTurbine%CtCurve(2, i)
!		print*, inTurbine%CtCurve(1, i), inTurbine%CpCurve(2, i),inTurbine%CtCurve(2, i)
!	end do				
!	close(1001)
!	stop
				
	return
end subroutine readTurbineData

subroutine readTerrain(filename)
	
	implicit none
	character(*) :: filename
	integer i, j, k, nline, ierr
	real(8) :: bbox(2,3), lcorner(3), ucorner(3)
	real(8), allocatable :: z(:)
	real(8) :: tmpr, x, y, zr
	
	

	if( trim(terrainFileCurrent) == trim(filename)  ) then
		return;
	else
		if(allocated(terrain)) deallocate(terrain)
		if(allocated(octree_RBF)) deallocate(octree_RBF)
		call clean_tree	
		terrainFileCurrent = trim(filename)
	end if
!	if(allocated(terrain)) return



	open(unit = 101, file = trim(filename), status = 'unknown')
	nline = 0
	do i = 1, 100000
		read(101,*,iostat=ierr) tmpr
		if(ierr < 0) go to 100
		nline = nline + 1
	end do
100	continue

	lcorner = 10.0d10
	ucorner = 0.0d0
	rewind(101)
	if(.not. allocated(terrain)) allocate(terrain(nline,3))
	if(.not. allocated(z)) allocate(z(nline))
	do i = 1, nline
		read(101,*,iostat=ierr) terrain(i,1:2), z(i)
		terrain(i,3) = 0.0d0
		do j = 1, 3
			if(terrain(i,j) < lcorner(j)) lcorner(j) = terrain(i,j)
			if(terrain(i,j) > ucorner(j)) ucorner(j) = terrain(i,j)
		end do

		if(ierr < 0) go to 101
	end do
101	continue
	close(101)
	
	
	bbox = 0.0d0
	bbox(1,:) = lcorner - 10.0d0
	bbox(2,:) = ucorner	+ 10.0d0
	
	call octree_init(octreeMaxPts, octreeDepth, bbox)
	call octree_build([(i, i=1,nline)],terrain)
	call find_all_neighbors_recursive(tree%root_node)
	
	terrain(:,3) = z
	deallocate(z)

	if(trim(WindFLO_input%terrainModel) == 'RBF') then
		call constructOctreeRBF
	end if
	

	
!	open(unit = 101, file = 'testing.dat', status = 'unknown')
!	open(unit = 102, file = 'interpolated_res.dat', status = 'unknown')
!	do i = 1, 40000
!		read(101,*) x, y, zr
!		call getElevation(x,y,tmpr)
!		write(102,*) x, y, zr, tmpr
!	end do
!	close(101)
!	close(102)
!	stop
!	call get_all_depths(tree%root_node)	
!	call write_octree

	return
end subroutine readTerrain


subroutine readWindRose(filename)
	
	implicit none
	character(*) :: filename
	integer i, nline, ierr
	real(8) :: tmpr


	if( trim(windRoseFileCurrent) == trim(filename)  ) then
		return;
	else
		if(allocated(WindFLO_input%windRose)) deallocate(WindFLO_input%windRose)
		windRoseFileCurrent = trim(filename)
	end if	
!	if(allocated(WindFLO_input%windRose)) return
		
	open(unit = 101, file = trim(filename), status = 'unknown')
	nline = 0
	do i = 1, 100000
		read(101,*,iostat=ierr) tmpr
		if(ierr < 0) go to 100
		nline = nline + 1
	end do
100	continue
	rewind(101)

	if(.not. allocated(WindFLO_input%windRose)) allocate(WindFLO_input%windRose(nline))
	do i = 1, nline
		read(101,*) WindFLO_input%windRose(i)%modelVelocity, WindFLO_input%windRose(i)%probability
		WindFLO_input%windRose(i)%probability = WindFLO_input%windRose(i)%probability * 0.01d0
	end do
	close(101)
	
	return
end subroutine readWindRose



subroutine initializeTurbineProps(inTurbine)

	implicit none
	
	type(turbineProps) :: inTurbine

	inTurbine%turbineNum = undefined_i
	inTurbine%name = ''		
	inTurbine%position = undefined
	inTurbine%orientation = undefined
	inTurbine%height = undefined
	inTurbine%radius = undefined
	inTurbine%diameter = undefined
	inTurbine%ratedPower = undefined
	inTurbine%fictitious = .false.

	return
end subroutine initializeTurbineProps



subroutine Clean()  bind(C,name='Clean_')

	implicit none
	
	if(allocated(WindFLO_input%turbines)) deallocate(WindFLO_input%turbines)	
	if(allocated(turbines)) deallocate(turbines)
	
	return
end subroutine Clean

subroutine CleanAll()  bind(C,name='CleanAll_')

	implicit none
	
	call Clean
	if(allocated(WindFLO_input%windRose)) deallocate(WindFLO_input%windRose)
	if(allocated(terrain)) deallocate(terrain)
	if(allocated(octree_RBF)) deallocate(octree_RBF)
	call clean_tree
	
	return
end subroutine CleanAll


subroutine GetTurbine(ith, name, turbineParams , Ctx, Cty, Cpx, Cpy)  bind(C,name='GetTurbine_')
	
	implicit none

	integer ith	, i
	character(c_char) :: name(100)
	real(c_double) :: turbineParams(50)
	integer :: bool

	real(c_double) :: Ctx(5000), Cty(5000)
	real(c_double) :: Cpx(5000), Cpy(5000)
	
	
	turbineParams = 0.0d0;
	
	name = getCString(WindFLO_input%turbines(ith)%name)

	turbineParams(1) = WindFLO_input%turbines(ith)%turbineNum

	turbineParams(2) = WindFLO_input%turbines(ith)%position(1)
	turbineParams(3) = WindFLO_input%turbines(ith)%position(2)
	turbineParams(4) = WindFLO_input%turbines(ith)%position(3)
		
	turbineParams(5) = WindFLO_input%turbines(ith)%height
	turbineParams(6) = WindFLO_input%turbines(ith)%radius
	turbineParams(7) = WindFLO_input%turbines(ith)%diameter
	
	bool = WindFLO_input%turbines(ith)%fictitious
	turbineParams(8) = bool
	turbineParams(9) = WindFLO_input%turbines(ith)%ratedPower		
	
	turbineParams(10) = WindFLO_input%turbines(ith)%orientation(1)
	turbineParams(11) = WindFLO_input%turbines(ith)%orientation(2) 
	turbineParams(12) = WindFLO_input%turbines(ith)%orientation(3)
	
	bool = WindFLO_input%turbines(ith)%yaw	
	turbineParams(13) = bool
	
	do i = 1, size(WindFLO_input%turbines(ith)%CtCurve(1,:))
		Ctx(i) = WindFLO_input%turbines(ith)%CtCurve(1,i)
		Cty(i) = WindFLO_input%turbines(ith)%CtCurve(2,i)
	end do
	Ctx(i) = undefined
	Cty(i) = undefined	


	do i = 1, size(WindFLO_input%turbines(ith)%CpCurve(1,:))
		Cpx(i) = WindFLO_input%turbines(ith)%CpCurve(1,i)
		Cpy(i) = WindFLO_input%turbines(ith)%CpCurve(2,i)
	end do
	Cpx(i) = undefined
	Cpy(i) = undefined	
	
	
	return
end subroutine GetTurbine




subroutine GetRunParameters(runParams) bind(C,name='GetRunParameters_')

	implicit none
	real(c_double) :: runParams(50)
	integer bool

	runParams = 0.0d0;
	bool = WindFLO_input%batch
	runParams(1) = WindFLO_input%nTurbines
	runParams(2) = bool 
	if(allocated(WindFLO_input%windRose)) then
		runParams(3) = size(WindFLO_input%windRose)
	else
		runParams(3) = 0	
	end if

	return
end subroutine GetRunParameters

!subroutine GetAtmParameters(rho, turbulenceIntensity) bind(C,name='GetAtmParameters_')
subroutine GetAtmParameters(atmParams) bind(C,name='GetAtmParameters_')

	implicit none
	real(c_double) :: atmParams(50)

	atmParams = 0.0d0;
	atmParams(1) = WindFLO_input%rho
	atmParams(2) = WindFLO_input%turbulenceIntensity

	return
end subroutine GetAtmParameters


!subroutine GetAmbientWindParameters(modelName, modelVelocity,surfaceRoughness , referenceHeight)&
!&									bind(C,name='GetAmbientWindParameters_')
subroutine GetAmbientWindParameters(modelName, ambientParams) bind(C,name='GetAmbientWindParameters_')

	implicit none
	character(c_char) :: modelName(100)
	real(c_double) :: ambientParams(50)

	modelName = getCString(WindFLO_input%windModel)

	ambientParams = 0.0d0;
	ambientParams(1:3) = WindFLO_input%modelVelocity
	ambientParams(4) = WindFLO_input%surfaceRoughness
	ambientParams(5) = WindFLO_input%referenceHeight

	return
end subroutine GetAmbientWindParameters


!subroutine GetWakeParameters(wakeModelName, wakeMergeModelName , wakeExpansionCoeff) bind(C,name='GetWakeParameters_')
subroutine GetWakeParameters(wakeModelName, wakeMergeModelName , wakeParams) bind(C,name='GetWakeParameters_')

	implicit none
	character(c_char) :: wakeModelName(100)
	character(c_char) :: wakeMergeModelName(100)	
	real(c_double) :: wakeParams(50)

	wakeParams = 0.0d0;
	wakeModelName = getCString(trim(WindFLO_input%wakeModel))
	wakeMergeModelName = getCString(trim(WindFLO_input%wakeMergeModel))
	wakeParams(1:2) = WindFLO_input%wakeExpansionCoeff

	return
end subroutine GetWakeParameters


!subroutine GetNumericalParameters(gaussOrder, monteCarloPts) bind(C,name='GetNumericalParameters_')
subroutine GetNumericalParameters(numericalParams) bind(C,name='GetNumericalParameters_')

	implicit none
	integer(c_int) :: numericalParams(50)

	numericalParams = 0.0d0;
	numericalParams(1) = WindFLO_input%gaussOrder
	numericalParams(2) = WindFLO_input%monteCarloPts	

	return
end subroutine GetNumericalParameters


!subroutine GetNumericalParameters(gaussOrder, monteCarloPts) bind(C,name='GetNumericalParameters_')
subroutine GetCostParameters(costParams) bind(C,name='GetCostParameters_')

	implicit none
	real(c_double) :: costParams(50)

	costParams = 0.0d0;
	costParams(1) = WindFLO_input%coe
	
	return
end subroutine GetCostParameters



subroutine GetElevation(x,y,z)  bind(C,name='GetElevation_')

	implicit none
	real(c_double) :: x, y, z
			
	if(.not. allocated(terrain)) return		
			
	if(trim(WindFLO_input%terrainModel) == 'IDW') then
		call computeElevation_IDW((/x,y/),z)
	else if(trim(WindFLO_input%terrainModel) == 'RBF') then
		call computeElevation_RBF((/x,y/), z)
	else
		print*, 'Incorrect Terrain Model: Stopping'
		stop
	end if
	
	return
end subroutine GetElevation



subroutine GetWindRoseParameters(ith, windRoseParams) bind(C, name = 'GetWindRoseParameters_')

	implicit none
	integer(c_int) :: ith
	real(c_double) :: windRoseParams(50)
	real(8) :: v1(3), v2(3), dV, dTheta
	
	if(.not. allocated(WindFLO_input%windRose)) return

	windRoseParams(1:3) = WindFLO_input%windRose(ith)%modelVelocity
	windRoseParams(4) = WindFLO_input%windRose(ith)%probability	

	if(ith == 1) then
		v1 = WindFLO_input%windRose(ith)%modelVelocity
		v2 = WindFLO_input%windRose(ith+1)%modelVelocity
	else
		v1 = WindFLO_input%windRose(ith)%modelVelocity
		v2 = WindFLO_input%windRose(ith-1)%modelVelocity
	end if

	dV = abs(norm2(v1) - norm2(v2))
	dTheta =  abs(( 180.0d0 / pi) * acos(dot_product(v1,v2) / (norm2(v1)*norm2(v2))))

	windRoseParams(5) = 1.0d0 ! dv*dTheta


	return
end subroutine GetWindRoseParameters



subroutine ComputeRotationMatrix(a, b, R) bind(C,name='GetRotationMatrix_')

	implicit none
	real(8) :: a(3), b(3), R(3,3)
	real(8) :: u(3,1)
	real(8) :: x(3), v(3,3)
	real(8) :: theta, eye(3,3)
	integer :: i(1)
	
	eye = 0.0d0;
	eye(1,1) = 1.0d0; eye(2,2) = 1.0d0; eye(3,3) = 1.0d0
	
	x = cross(a,b)
	x = x/norm2(x)
	theta = acos ( (dot_product(a,b)) / (norm2(a) * norm2(b)))


	if(theta < 1.0d-6) then
		R = 0.0d0
		R(1,1) = 1.0d0
		R(2,2) = 1.0d0
		R(3,3) = 1.0d0
		
		WindFLO_input%rotMat = R
		WindFLO_input%rotMat_inv = R
		call MatInv3(WindFLO_input%rotMat_inv)
		
		return
	else if( (pi - theta) < 1.0d-2) then
		i = minloc(abs(a))
		x = 0.d0; x(i(1)) = 1.0d0
		x = cross(a,x)
		x = x/norm2(x)
	end if

	
	u(:,1) = x
	
	v(1,1) = 0.0d0; v(1,2) = -x(3) ; v(1,3) = x(2)
	v(2,1) = x(3); v(2,2) = 0.0d0 ; v(2,3) = -x(1)
	v(3,1) = -x(2); v(3,2) = x(1) ; v(3,3) = 0.0d0
		
	R = eye * cos(theta) + sin(theta)*v + (1.0d0 - cos(theta))*matmul(u , transpose(u))
		
	WindFLO_input%rotMat = R
	WindFLO_input%rotMat_inv = R
	
	call MatInv3(WindFLO_input%rotMat_inv)
				
	return
end subroutine ComputeRotationMatrix


subroutine ComputeArcLength(A, B, length) bind(C,name='ComputeArcLength_')

	implicit none
	integer, parameter :: N = 500
	real(8) :: A(3), B(3), length
	real(8) :: old(3), new(3), stepV(3), dist
	integer i
	


	stepV = 0.0d0; stepV(1) = (B(1) - A(1)) / (1.0d0 * N)
	
	
	old = matmul(WindFLO_input%rotMat_inv, A)
	call GetElevation(old(1),old(2),old(3))
	dist = 0.0d0
	length = 0.0d0
	do i = 1, N
		new = A + stepV * i
		new = matmul(WindFLO_input%rotMat_inv, new)
		call GetElevation(new(1),new(2),new(3))
		if(new(3) /= undefined) then
			dist = norm2(new-old)
			old = new		
		end if
		
		length = length + dist
	end do





!	old = matmul(WindFLO_input%rotMat_inv, A)
!	new = matmul(WindFLO_input%rotMat_inv, B)
!	call GetElevation(old(1),old(2),old(3))
!	call GetElevation(new(1),new(2),new(3))
!	stepV = (new - old)/ (1.0d0 * N)
!	do i = 1, N
!		new = old + stepV
!		call GetElevation(new(1),new(2),new(3))
!		length = length + norm2(new-old)
!		old = new		
!	end do


!	do i = 1, N
!		new = A
!		new(1) = A(1) + step * i		
!		new = matmul(WindFLO_input%rotMat_inv, new)
!		call GetElevation(new(1),new(2),new(3))
!		length = length + norm2(new-old)
!		old = new		
!	end do
	
	return
end subroutine ComputeArcLength







character(len=1000) function getFString(inString)
	implicit none
	character(c_char) :: inString(:)
	integer i, n
			
	n = size(inString)		
	do i = 1, n
		if(inString(i) == '\0' .or. inString(i) == '' .or. inString(i) == ' ') then
			go to 100
		else
			getFString(i:i) = inString(i)
		end if				
	end do
100	continue	

	return
end function getFString


function getCString(inString) result(outString)
	implicit none
	character(*) :: inString	
	character(c_char) :: outString(100)
	integer i, n
			
	n = len(inString)
	
	do i = 1, n
	outString(i) = ''
		if(inString(i:i) /= '\0' .and. inString(i:i) /= '') then
			outString(i) = inString(i:i)
		else
			goto 100
		end if				
	end do
100	continue	

	outString(i) = C_NULL_CHAR
	
	return
end function getCString



subroutine getCpOrCtCurve(CpCurve,CtCurve)

	implicit none
	
	real(8) :: CpCurve(2,5000)			! Cp curve data for spline fitting
	real(8) :: CtCurve(2,5000)			! Ct curve data for spline fitting
	real(8) :: a, coeff(4), x(3)
	
	integer i


	do i = 1, 5000
	

		if( CtCurve(1,i) /= undefined .and. CtCurve(2,i) /= undefined .and. &
		&   CpCurve(1,i) == undefined .and. CpCurve(2,i) == undefined ) then
	
			coeff(1) = 0.0d0
			coeff(2) = -4.0d0
			coeff(3) = 4.0d0
			coeff(4) = -CtCurve(2,i)

			call ccubsolv(coeff,x)
			a = minval(x)

			CpCurve(1,i) = CtCurve(1,i)
			CpCurve(2,i) = (4.0d0 * a*(1.0d0 - a)**2.0d0)
		
		else if( CtCurve(1,i) == undefined .and. CtCurve(2,i) == undefined .and. &
		&   	 CpCurve(1,i) /= undefined .and. CpCurve(2,i) /= undefined ) then
				
			coeff(1) = 4.0d0
			coeff(2) = -8.0d0
			coeff(3) = 4.0d0
			coeff(4) = -CpCurve(2,i)
		
			call ccubsolv(coeff,x)
			a = minval(x)
			
			CtCurve(1,i) = CpCurve(1,i)
			CtCurve(2,i) = (4.0d0 * a*(1.0d0 - a))

		else		
			return
		end if		
	 	
	end do

	return		
end subroutine getCpOrCtCurve





function cross(a,b) result(c)
	implicit none
	
	real(8) :: a(3), b(3), c(3)
	c = (/a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3) , a(1)*b(2) - a(2)*b(1) /) 

	return
end function cross

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+
!
!     coeff(1)*x^3 + coeff(2)*x^2 + coeff(3)*x + coeff(4) = 0
!	  x^3 + a*x^2 + b*x + c = 0
!
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+
	subroutine ccubsolv(coeff,x)
	real(8) :: coeff(4),x(3)
	real(8) :: pi
	real(8) :: a, b, c, Q, R
	real(8) :: theta
	
	pi=3.14159265358979323846d0
 
 	if(coeff(1) /= 0.0d0) then
 
	 	a = coeff(2) / coeff(1)
	 	b = coeff(3) / coeff(1)
	 	c = coeff(4) / coeff(1)
 
	 	Q = ( a*a - 3.0d0*b ) / (9.0d0)
 		R = (2.0d0*a*a*a - 9.0d0 *a*b + 27.0d0*c) / (54.0d0)
 
	 	theta = acos(R / sqrt(Q**3.0d0))
 	
	 	x(1) = -2.0d0 * sqrt(Q) * cos(theta/ 3.0d0) - (a / 3.0d0)
	 	x(2) = -2.0d0 * sqrt(Q) * cos((theta + 2.0d0*pi)/ 3.0d0) - (a / 3.0d0)
	 	x(3) = -2.0d0 * sqrt(Q) * cos((theta - 2.0d0*pi)/ 3.0d0) - (a / 3.0d0)
 
 	else
	 	a = coeff(2) 
	 	b = coeff(3) 
	 	c = coeff(4) 
 		
 		Q = -0.5d0 * (b + sign(-1.0d0,b) * sqrt(b*b - 4.0d0*a*c) )
	 	x(1) = Q / a
	 	x(2) = c / Q
 		x(3) = undefined
 	
 	end if
 
	return
	
	end subroutine ccubsolv



subroutine MatInv3(A) bind(C,name='MatInv3_')

	implicit none

    real(8) :: A(3,3)
    real(8) :: B(3,3)
    real(8) :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    
    A = B
    
    return
  end subroutine MatInv3


end module WindFLO



