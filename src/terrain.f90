!
!
!	Copyright 2019, Sohail R. Reddy
!	email: sredd001@fiu.edu
!	www.sohailreddy.com
!
!

module terrain_mod

	use iso_c_binding
	use octree_mod
	
	type RBF
		logical :: constructed
		integer, allocatable :: ids(:)
		real(8), allocatable :: w(:)
	end type RBF
	
	
	type(RBF), allocatable :: octree_RBF(:)
	
	integer :: rbfKernel = 1
	integer :: nNearest = 10	! NOT USED ANYMORE
	real(8) :: shapeFactor = 1.0d0
	real(8) :: powerIDW = 3.0d0
	integer :: octreeDepth = 10
	integer :: octreeMaxPts = 10
	
	contains



! Inverse Distance Weighting
subroutine computeElevation_IDW(x_new,z)

	implicit none
	real(8) :: x_new(:), z
	
	integer :: i, j, k
	real(8) :: x, y, dist, power
	type(node_type), pointer :: cnode
	integer :: ids(tree%max_pts*10)
	real(8) :: weights(tree%max_pts*10)
	real(8), parameter :: undefined = 9876543210.0d0	


	x = x_new(1)
	y = x_new(2)
			
	power = powerIDW
	nullify(cnode)
	call octree_search((/x,y,0.0d0/), 0.0d0, cnode = cnode)

	if(.not. associated(cnode)) then
		z = undefined
		return
	end if


	weights = 0.0d0
	
	j = 0
	ids = 0
	if(cnode%num_point > 0) then
		do i = 1, cnode%num_point
			j = j + 1
			ids(j) = cnode%point_ids(i)
		end do
	end if
	
	if(associated(cnode%neighbors)) then
		do i = 1, size(cnode%neighbors)
			do k = 1, cnode%neighbors(i)%num_point
				j = j + 1
				ids(j) = cnode%neighbors(i)%point_ids(k)
!				if( j >= nNear) go to 100
			end do
		end do
	end if	
100	continue		

	do i = 1, j
		dist = norm2( tree%points(ids(i),1:2) - (/x,y/) )
		if(dist == 0.0d0) then
			z = tree%points(ids(i),3)
			return
		end if
		weights(i) = 1.0d0 / (dist ** power)
	end do

	weights(1:j) = weights(1:j)/ sum(weights(1:j))
	z = dot_product(weights(1:j),  tree%points(ids(1:j),3))

	return
end subroutine computeElevation_IDW



!	Radial Basis Functions
subroutine computeElevation_RBF(x_new, z)

	implicit none
	real(8) :: x_new(:), z
	
	integer :: id
	real(8) :: x, y
	type(node_type), pointer :: cnode
	real(8) :: c = 1.0d0
	real(8), parameter :: undefined = 9876543210.0d0	


	x = x_new(1)
	y = x_new(2)
	
	nullify(cnode)			
	call octree_search((/x,y,0.0d0/), 0.0d0, cnode = cnode)

	if(.not. associated(cnode)) then
		z = undefined
		return
	end if

	id = cnode%leaf_id 
	if(	octree_RBF(id)%constructed) then
		z = evaluateRBF(tree%points( octree_RBF(id)%ids , 1:2), x_new, octree_RBF(id)%w, c)
	else
		print*, 'RBF Not Constructed: Using IDW'
		call computeElevation_IDW(x_new,z)
!		stop
	end if
	return
end subroutine computeElevation_RBF



subroutine constructOctreeRBF

	implicit none

	integer :: i, j, k, npts
	type(node_type), pointer :: cnode
	
	integer :: ids(tree%max_pts*9)
	real(8) :: c = 1.0d0


	call get_all_leaves
	if(.not. allocated(octree_RBF)) allocate(octree_RBF(tree%nleaves))
	

	do i = 1, tree%nleaves
	
		cnode => tree%leaves(i)%node
		cnode%leaf_id = i
	
		npts = 0
		ids = 0

		if(cnode%num_point > 0) then
			do j = 1, cnode%num_point
				npts = npts + 1
				ids(npts) = cnode%point_ids(j)
			end do
		end if
	
		if(associated(cnode%neighbors)) then
			do j = 1, size(cnode%neighbors)
				if(cnode%neighbors(j)%num_point > 0) then
					do k = 1, cnode%neighbors(j)%num_point
						npts = npts + 1
						ids(npts) = cnode%neighbors(j)%point_ids(k)
					end do
				end if
			end do
		end if	

		if(npts == 0) cycle

		if(.not. allocated(octree_RBF(i)%ids)) allocate(octree_RBF(i)%ids(npts))
		if(.not. allocated(octree_RBF(i)%w)) allocate(octree_RBF(i)%w(npts))
		
		octree_RBF(i)%ids = ids(1:npts)
		call constructRBF(tree%points(ids(1:npts),1:2), tree%points(ids(1:npts),3) ,octree_RBF(i)%w,c)
		octree_RBF(i)%constructed = .true.

	end do

	return
end subroutine constructOctreeRBF


subroutine constructRBF(xyz, f, w,c)

	implicit none
	real(8) :: xyz(:,:)
	real(8) :: f(:), w(:)
	real(8), allocatable :: k(:,:)
	real(8) :: c 
	
	integer i, j, n
	
	w = 0.0d0
	
	n = size(w)
	allocate(k(n,n))
	
	do i = 1, n
		do j = i , n
			k(i,j) = rbfKernelEval(norm2(xyz(i,:) - xyz(j,:)))
			k(j,i) = k(i,j)
		end do
	end do
	
	call  solve_lin_eq(n, k, f, w)
		
	return
end subroutine constructRBF

function evaluateRBF(xyz, x, w, c) result(f)

	implicit none
	real(8) :: xyz(:,:)
	real(8) :: w(:), x(:)
	real(8) :: f, c
	integer i, j, n
	
	f = 0.0d0
	n = size(w)
	do i = 1, n
		f = f + w(i) * rbfKernelEval(norm2(xyz(i,:) - x))
	end do

	return
end function evaluateRBF



function rbfKernelEval(d) result(ans)

	implicit none
	real(8) :: d, ans
	real(8) :: c
	c = shapeFactor
	if(rbfKernel == 1) then
		ans = sqrt(d*d + c*c)	
	else if (rbfKernel == 2) then
		ans = 1.0d0/sqrt(d*d + c*c)	
	else if (rbfKernel == 3) then
		ans = exp(-(d*d) / (c*c))
	end if	

	return
end function rbfKernelEval

end module terrain_mod






subroutine solve_lin_eq(neq, a, b, x)

	implicit none      
	double precision       a(neq, neq)
	double precision       b(neq)
	double precision       x(neq)
	integer                indx(neq)
	double precision       dd
	integer                i, neq

	call svd2(neq,a,x,b)
	
	do i = 1, neq
		if (x(i) /= x(i) ) x(i) = 1.0d0
	enddo

	return
end subroutine solve_lin_eq





module svdmodule
  implicit none
  integer mp,np,m,n
  real(8),dimension(:,:),allocatable :: v,u
  real(8),dimension(:),allocatable   :: w,c
end module svdmodule

subroutine svd2(n1,a,x,b)
  use svdmodule
  implicit none
  integer k,n1
  real(8) a(n1,n1),x(n1),b(n1)
  real(8) wmax,wmin,rms
!
  mp=n1
  np=n1
  m=n1
  n=n1
  allocate(v(n,n),u(n,n))
  allocate(w(n),c(n))
! copy a into u
  u=a
! decompose matrix a
  call svdcmp
! find maximum singular value
  wmax=0.d0
  do k=1,n
    if(w(k).gt.wmax)wmax=w(k)
  end do
! define "small"
  wmin=1.d-6 !wmax*tiny(wmax) !(1.d-6)
! zero the "small" singular values
  do k=1,n
    if(w(k).lt.wmin)w(k)=0.d0
  end do
! backsubstitute for each right-hand side vector
  c=b
  call svbksb(x)
  !rms=dsqrt(sum((matmul(a,x)-b)**2))
  deallocate(v,u)
  deallocate(w,c)
end subroutine svd2

subroutine svbksb(x)
  use svdmodule
  implicit none
  real(8) x(np)
  integer i,j,jj
  real(8) s,tmp(n)
  do j=1,n
    s=0.d0
    if(w(j).ne.0.d0)then
      do i=1,m
        s=s+u(i,j)*c(i)
      end do
      s=s/w(j)
    endif
    tmp(j)=s
  end do
  do j=1,n
    s=0.d0
    do jj=1,n
      s=s+v(j,jj)*tmp(jj)
    end do
    x(j)=s
  end do
end subroutine svbksb

subroutine svdcmp
  use svdmodule
  implicit none
  integer i,its,j,jj,k,l,nm
  real(8) anorm,c1,f,g,h,s,scale,x,y,z,rv1(n),pythag2
  g=0.d0
  scale=0.d0
  anorm=0.d0
  do i=1,n
    l=i+1
    rv1(i)=scale*g
    g=0.d0
    s=0.d0
    scale=0.d0
    if(i.le.m)then
      do k=i,m
        scale=scale+dabs(u(k,i))
      end do
      if(scale.ne.0.d0)then
        do k=i,m
          u(k,i)=u(k,i)/scale
          s=s+u(k,i)*u(k,i)
        end do
        f=u(i,i)
        g=-dsign(dsqrt(s),f)
        h=f*g-s
        u(i,i)=f-g
        do j=l,n
          s=0.d0
          do k=i,m
            s=s+u(k,i)*u(k,j)
          end do
          f=s/h
          do k=i,m
            u(k,j)=u(k,j)+f*u(k,i)
          end do
        end do
        do k=i,m
          u(k,i)=scale*u(k,i)
        end do
      endif
    endif
    w(i)=scale *g
    g=0.d0
    s=0.d0
    scale=0.d0
    if((i.le.m).and.(i.ne.n))then
      do k=l,n
        scale=scale+dabs(u(i,k))
      end do
      if(scale.ne.0.d0)then
        do k=l,n
          u(i,k)=u(i,k)/scale
          s=s+u(i,k)*u(i,k)
        end do
        f=u(i,l)
        g=-dsign(dsqrt(s),f)
        h=f*g-s
        u(i,l)=f-g
        do k=l,n
          rv1(k)=u(i,k)/h
        end do
        do j=l,m
          s=0.d0
          do k=l,n
            s=s+u(j,k)*u(i,k)
          end do
          do k=l,n
            u(j,k)=u(j,k)+s*rv1(k)
          end do
        end do
        do k=l,n
          u(i,k)=scale*u(i,k)
        end do
      endif
    endif
    anorm=max(anorm,(dabs(w(i))+dabs(rv1(i))))
  end do
  do i=n,1,-1
    if(i.lt.n)then
      if(g.ne.0.d0)then
        do j=l,n
          v(j,i)=(u(i,j)/u(i,l))/g
        end do
        do j=l,n
          s=0.d0
          do k=l,n
            s=s+u(i,k)*v(k,j)
          end do
          do k=l,n
            v(k,j)=v(k,j)+s*v(k,i)
          end do
        end do
      endif
      do j=l,n
        v(i,j)=0.d0
        v(j,i)=0.d0
      end do
    endif
    v(i,i)=1.d0
    g=rv1(i)
    l=i
  end do
  do i=min(m,n),1,-1
    l=i+1
    g=w(i)
    do j=l,n
      u(i,j)=0.d0
    end do
    if(g.ne.0.d0)then
      g=1.d0/g
      do j=l,n
        s=0.d0
        do k=l,m
          s=s+u(k,i)*u(k,j)
        end do
        f=(s/u(i,i))*g
        do k=i,m
          u(k,j)=u(k,j)+f*u(k,i)
        end do
      end do
      do j=i,m
        u(j,i)=u(j,i)*g
      end do
    else
      do j= i,m
        u(j,i)=0.d0
      end do
    endif
    u(i,i)=u(i,i)+1.d0
  end do
  do k=n,1,-1
    do its=1,30
      do l=k,1,-1
        nm=l-1
        if((dabs(rv1(l))+anorm).eq.anorm)  goto 2
        if((dabs(w(nm))+anorm).eq.anorm)  goto 1
      end do
1     c1=0.d0
      s=1.d0
      do i=l,k
        f=s*rv1(i)
        rv1(i)=c1*rv1(i)
        if((dabs(f)+anorm).eq.anorm) goto 2
        g=w(i)
        h=pythag2(f,g)
        w(i)=h
        h=1.d0/h
        c1= (g*h)
        s =-(f*h)
        do j=1,m
          y=u(j,nm)
          z=u(j,i)
          u(j,nm)=(y*c1)+(z*s)
          u(j,i)=-(y*s)+(z*c1)
        end do
      end do
2     z=w(k)
      if(l.eq.k)then
        if(z.lt.0.d0)then
          w(k)=-z
          do j=1,n
            v(j,k)=-v(j,k)
          end do
        endif
        goto 3
      endif
    !  if(its.eq.100) pause 'no convergence in svdcmp'  !30
      x=w(l)
      nm=k-1
      y=w(nm)
      g=rv1(nm)
      h=rv1(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
      g=pythag2(f,1.d0)
      f=((x-z)*(x+z)+h*((y/(f+dsign(g,f)))-h))/x
      c1=1.d0
      s=1.d0
      do j=l,nm
        i=j+1
        g=rv1(i)
        y=w(i)
        h=s*g
        g=c1*g
        z=pythag2(f,h)
        rv1(j)=z
        c1=f/z
        s =h/z
        f= (x*c1)+(g*s)
        g=-(x*s)+(g*c1)
        h=y*s
        y=y*c1
        do jj=1,n
          x=v(jj,j)
          z=v(jj,i)
          v(jj,j)= (x*c1)+(z*s)
          v(jj,i)=-(x*s)+(z*c1)
        end do
        z=pythag2(f,h)
        w(j)=z
        if(z.ne.0.d0)then
          z=1.d0/z
          c1=f*z
          s =h*z
        endif
        f= (c1*g)+(s*y)
        x=-(s*g)+(c1*y)
        do jj=1,m
          y=u(jj,j)
          z=u(jj,i)
          u(jj,j)= (y*c1)+(z*s)
          u(jj,i)=-(y*s)+(z*c1)
        end do
      end do
      rv1(l)=0.d0
      rv1(k)=f
      w(k)=x
    end do
3   continue
  end do
end subroutine svdcmp

real(8) function pythag2(a,b)
  implicit none
  real(8) a,b
  real(8) absa,absb
  absa=dabs(a)
  absb=dabs(b)
  if(absa.gt.absb)then
    pythag2=absa*dsqrt(1.d0+(absb/absa)**2)
  else
    if(absb.eq.0.d0)then
      pythag2=0.d0
    else
      pythag2=absb*dsqrt(1.d0+(absa/absb)**2)
    endif
  endif
end

