### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 398fdf3a-6854-11eb-3041-b9130d132483
using Random, Distributions

# ╔═╡ 3a5fc9dc-68ab-11eb-3a40-9fc97f9b2ebe
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add(["PlutoUI", "Plots"])
	using PlutoUI
end

# ╔═╡ ec133ca2-6763-11eb-002e-45a5d78c3605
md"
The First step is to gather the data from appendix and create a csv file out of it
"

# ╔═╡ c9c34aa6-675f-11eb-1fc5-196e34a58196
species = [ "Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria gallica",
			"Armillaria sinapina",
			"Armillaria tabescens",
			"Armillaria tabescens",
			"Fomes fomentarius",
			"Hypohodontia crustosa",
			"Hyphoderma setigerum",
			"Hyphoderma setigerum",
			"Laetiporus conifericola",
			"Lentinus crinitus",
			"Mycoacia meridionalis",
			"Merulius tremellosus",
			"Merulius tremellosus",
			"Phlebiopsis flavidoalba",
			"Phlebiopsis flavidoalba",
			"Phellinus glivus",
			"Phellinus hartigii",
			"Porodisculus pendulus",
			"Phellinus robiniaee",
			"Phellinus robiniaee",
			"Phlebia acerina",
			"Phlebia acerina",
			"Pycnoporus sanguineus",
			"Schizophyllum commune",
			"Schizophyllum commune",
			"Tyromyces chioneus",
			"Xylobolus subpileatus"]


# ╔═╡ b6da61c0-6775-11eb-1352-cf2f82d47fea
size(species)

# ╔═╡ 0e794288-676d-11eb-0d34-1fc7c09880d5
unique_species = unique(species)

# ╔═╡ 49420896-6761-11eb-025d-c7b05f2e0f08
DR = [  8.93,
	    4.54,
		5.81,
		6.95,
		5.03,
		2.04,
		3.44,
		2.00,
		4.58,
		1.89,
		5.10,
		21.87,
		9.21,
		7.61,
		6.57,
		7.07,
		8.02,
		4.60,
		34.01,
		21.15,
		17.95,
		11.38,
		16.16,
		8.01,
		3.02,
		5.17,
		7.10,
		16.78,
		20.98,
		14.49,
		4.65,
		3.88,
		13.75,
		5.17 ]

# ╔═╡ f6f8dca2-6764-11eb-1ea4-11949b018219
md" I calculate the geometric mean based on $(temp10 * temp16 * temp22)^{1/3}$"

# ╔═╡ 1dae4998-6765-11eb-3906-5158627996f5
HR_data = [	[0.30, 0.36, 0.34],
		[0.18, 0.26, 0.38],
		[0.26, 0.24, 0.32],
		[0.16, 0.30, 0.24],
		[0.20, 0.24, 0.40],
		[0.14, 0.32, 0.48],
		[0.20, 0.26, 0.36],
		[0.06, 0.18,  0.66],
		[0.33, 0.60, 0.84],
		[0.35, 0.60, 0.93],
		[0.32, 0.68, 1.56],
		[0.36, 1.28, 4.62],
		[1.20, 0.99, 1.77],
		[1.39, 3.70, 6.46],
		[0.44, 1.90, 3.68],
		[1.08, 3.31, 6.00],
		[1.64, 3.06, 6.17],
		[0.36, 1.10, 1.60],
		[3.30, 5.85, 8.67],
		[3.40, 6.50, 8.33],
		[2.28, 5.70, 8.41],
		[3.04, 7.40, 10.57],
		[1.40, 1.53, 3.70],
		[0.49, 1.26, 0.94],
		[0.95, 1.25, 2.90],
		[0.40, 1.52, 3.32],
		[0.39, 1.24, 2.84],
		[3.70, 7.40, 8.27],
		[3.70, 7.40, 8.23],
		[0.81, 3.21, 7.26],
		[1.88, 3.32, 7.40],
		[1.06, 1.64, 4.60],
		[1.92, 3.37, 5.67],
		[0.74, 1.00, 1.04]]

# ╔═╡ c4826416-6767-11eb-0e61-e912d529bebc
md" Clear up why are we using geometric mean  \
	Generate a table of our data
"

# ╔═╡ 54503b68-6768-11eb-369c-9d1714c8a21c
md" Now let's plot stuff "

# ╔═╡ 67f36fa8-676d-11eb-3055-bd1f93a5c936
scatter(log.(DR))

# ╔═╡ 7a738940-6776-11eb-39e7-c33e9301faf3
normalize(x) = (x.-minimum(x))/(maximum(x)-minimum(x))

# ╔═╡ 4c2ed13a-6851-11eb-3fcb-c79d02727c89
md"
# Simulating growth
"

# ╔═╡ 78c91270-6827-11eb-250a-6f5be3f36820
mutable struct Point
	x::Float64
	y::Float64
end

# ╔═╡ ec175d6c-6828-11eb-2f6f-c75c5b18c3ef
begin
	import Base.+
	(+)(p1::Point, p2::Point) = Point(p1.x+p2.x, p1.y+p2.y)
end

# ╔═╡ 6892328c-68ff-11eb-2a85-4f5aa9d53b32
begin
	import Base.*
	(*)(c::Float64, p::Point) = Point(p.x*c, p.y*c)
end

# ╔═╡ 62ff5ea8-6765-11eb-1ae9-e597834764df
geometric_mean(x) = (x[1]*x[2]*x[3])^(1/3)

# ╔═╡ a2efe0d2-6765-11eb-2249-c50e7ca75cd1
HR = geometric_mean.(HR_data)

# ╔═╡ 19726662-6761-11eb-1d7a-df42ac99ef1c
begin
	using DataFrames
	using CSV
	using Plots
	df = DataFrame(Species = species,
				   decomposition_rate = DR,
				   growth_rate = HR)
	#CSV.write("/data/MCM/dc_data.csv",df)
end

# ╔═╡ 452e74e2-6772-11eb-347f-3556436653ca
df2 = CSV.File("Fungal_trait_data.csv") |> DataFrame

# ╔═╡ a1719b38-6774-11eb-01fd-511391d46d16
begin
	data = df2["water.niche.width"][1:34]
	moisture_niche_width= normalize(data)
end

# ╔═╡ aadd29bc-6774-11eb-02c6-59f021e86c6a
competitive_ranking = normalize(df2["ranking"][1:34])

# ╔═╡ 35dfd5e6-6775-11eb-23d6-052efdcd1dd1
moisture_tolerance = (competitive_ranking - moisture_niche_width)

# ╔═╡ c7ecc5ca-6775-11eb-0da6-09045374100c
scatter(moisture_tolerance, log.(DR), legend=false)

# ╔═╡ 70cbb1c8-6768-11eb-36c4-b7f1e44ceb9d
scatter(HR, DR, 
		xaxis = "Hyphal Extension Rate", 
		yaxis="Decomposition Rate", 
		legend = false)

# ╔═╡ 8cea6466-6836-11eb-2704-797ca64396a1
begin
	df3 = DataFrame(Species = species,
				   decomposition_rate = DR,
				   growth_rate = HR,
				   moisture_tolerance = moisture_tolerance)
	CSV.write("/data/MCM/dc_data.csv",df)
end

# ╔═╡ 7af43442-677a-11eb-2a4f-71ca69516b35
begin
	df_final = DataFrame(Species = species,
				   decomposition_rate = DR,
				   growth_rate = HR, 
				   moisture_tolerance = moisture_tolerance)

	CSV.write("/data/MCM/dc_data2.csv",df_final)
end

# ╔═╡ e4efbe66-677a-11eb-105b-819ebbd28990
scatter(HR, DR, markersize=(moisture_tolerance.+1)*5, alpha=0.5,
		xaxis = "Hyphal Extension Rate", 
		yaxis="Decomposition Rate", 
		legend = false)

# ╔═╡ ed19a464-68ea-11eb-1d13-cd2480f4186c
function isequal(p1::Point,p2::Point)
	return p1.x== p2.x && p1.y==p2.y
end

# ╔═╡ 1ea0c1a8-683d-11eb-0efa-89c527035050
struct Hypha
	points::Vector{Point}
	isActive::Bool
end

# ╔═╡ 1c533af0-6829-11eb-2f2d-ab894ae524ca
filter(row -> row["Species"] == "Fomes fomentarius", df3)["decomposition_rate"][1]

# ╔═╡ 4d04ef40-68b0-11eb-2781-c9d566b2ece5
begin a = [Point(1,1)]
	append!(a,[Point(1,1)])
end

# ╔═╡ 55f94e50-6826-11eb-2a77-713890850b60
begin
mutable struct Fungi
	species::String
	center::Point
	death_rate::Float64
	speed::Float64
	extension_rate::Float64
	decomposition_rate::Float64
	moisture_tolerance::Float64
	from::Vector{Point}
	to::Vector{Point}
	angle::Vector{Float64}
	isTip::Vector{Bool}
	Si::Vector{Float64}


end

	Fungi(species::String,
		  center::Point,
		  death_rate::Float64,
		  speed::Float64) = Fungi(
					species,
					Point(1,1),
				    1.0,
				    1.0,
					filter(row -> row["Species"] .== species, df3)["decomposition_rate"][1],
					filter(row -> row["Species"] .== species, df3)["growth_rate"][1],
					filter(row -> row["Species"] .== species, df3)["moisture_tolerance"][1],
					[center,center,center,center,center,center,center,center],
					[1e-4*Point(1,0)+center,
					 1e-4*Point(sqrt(2)/2,sqrt(2)/2)+center,
					 1e-4*Point(-1,0)+center,
					 1e-4*Point(-sqrt(2)/2,sqrt(2)/2)+center,
					 1e-4*Point(0,1)+center,
					 1e-4*Point(-sqrt(2)/2,-sqrt(2)/2)+center,
					 1e-4*Point(0,-1)+center,
					 1e-4*Point(sqrt(2)/2,-sqrt(2)/2)+center,],
					[0,
					 pi/4,
					 pi/2,
					 3*pi/4,
					 pi,
					 pi+pi/4,
					 pi+pi/2,
					 pi+3*pi/4],
					 [true,
					  true,
					  true,
					  true,
					  true,
					  true,
					  true,
					  true],
					 [ 10^(-5)
					   10^(-5)
			    	   10^(-5)
					   10^(-5)
					   10^(-5)
					   10^(-5)
					   10^(-5)
					   10^(-5)]
					 )
end

# ╔═╡ 2674f916-68d2-11eb-0066-1fcc4722e530


# ╔═╡ d1293282-682f-11eb-1e05-71185e0dcc8b
 genk = Fungi("Fomes fomentarius",
		Point(50,50),
		1.0,
		1.0)

# ╔═╡ df20ac7a-6830-11eb-077f-5d0b553551c2
function getx(p::Point)
	return p.x
end

# ╔═╡ f6aea732-6830-11eb-1fbe-bfd1f16815a3
function gety(p::Point)
	return p.y
end

# ╔═╡ 01e89ad2-684e-11eb-1a44-fde170bc69d3
function bernoulli(p)
	return rand()<p
end

# ╔═╡ 8babf7d4-684e-11eb-1cc5-8d13e27b32b3
md" 
#   JUSTIFYYYYYYYYYYY 
"

# ╔═╡ b9127852-68f5-11eb-0b11-172f8bb57677
t = 0.01 # days

# ╔═╡ 7cc884ac-68fc-11eb-0500-b3e0622146dc
Se = zeros(10000,10000)

# ╔═╡ 469a53c2-68fe-11eb-1064-8d806c382970


# ╔═╡ 7651344a-68f5-11eb-17be-d7573461d1a8
function growth_rate()
	return 0.5
end

# ╔═╡ ac8dd314-6847-11eb-3f2e-1163d16252d6
function next_point(p::Point, ϕ::Float64, isTip::Bool, Si::Float64)
	
	if(!isTip)
		return -1
	end
	
	if(Si < 10^(-12)) 
		return 1
	end
	

	if(bernoulli(1/8))
		if(bernoulli(0.50))
			ϕ += π/6
		else
			ϕ -= π/6	
		end
	end
	v = growth_rate() 
	return (p, Point(p.x + v*t*cos(ϕ), p.y +  v*t*sin(ϕ)),ϕ)
		
end

# ╔═╡ e9e11d3c-6853-11eb-0b2f-c310a84e7f9e
function next_branch(p::Point)
	d=Normal(1.4,0.17)
	ϕ = rand(d)
	return Point(p.x + cos(ϕ), p.y +  sin(ϕ))
end

# ╔═╡ f3df1314-68b4-11eb-2d92-6b515469e676
function next_step(p::Point)
	return Point(p.x+1,p.y+1)
end

# ╔═╡ d5796fb4-68b9-11eb-021a-e537d93880d7
function getFirst(t::Tuple{Point,Point,Float64})
	return t[1]
end

# ╔═╡ f6034e94-68b9-11eb-3586-d18d276bd9be
function getSecond(t::Tuple{Point,Point,Float64})
	return t[2]
end

# ╔═╡ 331a2554-68bb-11eb-0fb9-99fc2c426c1a
function getThird(t::Tuple{Point,Point,Float64})
	return t[3]
end

# ╔═╡ c097e286-68b9-11eb-2799-333f9b0cd757
genk

# ╔═╡ d6eac3f2-68be-11eb-026d-6729ab9283b1
s = next_point.(genk.to,genk.angle,genk.isTip,genk.Si)


# ╔═╡ 1c596458-68c8-11eb-3993-29610d466d9a
function orientation(p, q, r)
    # to find the orientation of an ordered triplet (p,q,r) 
    # function returns the following values: 
    # 0 : Colinear points 
    # 1 : Clockwise points 
    # 2 : Counterclockwise 
      
    # See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/  
    # for details of below formula.  
      
    val = ((q.y - p.y) * (r.x - q.x)) - ((q.x - p.x) * (r.y - q.y)) 
    if (val > 0) 
        # Clockwise orientation 
        return 1
	elseif(val < 0) 
        return 2
    else 
        return 0
	end
end

# ╔═╡ 00acdf7c-68d3-11eb-273a-61da7ca2de6f
function onSegment(p, q, r) 
    if ( (q.x <= max(p.x, r.x)) && (q.x >= min(p.x, r.x)) && 
           (q.y <= max(p.y, r.y)) && (q.y >= min(p.y, r.y))) 
        return true
	end
    return false
end

# ╔═╡ 6d712a94-68c6-11eb-2881-1b394bde84a4
function intersects(p1, q1, p2, q2)     
    # Find the 4 orientations required for  
    # the general and special cases 
    o1 = orientation(p1, q1, p2) 
    o2 = orientation(p1, q1, q2) 
    o3 = orientation(p2, q2, p1) 
    o4 = orientation(p2, q2, q1) 
  	if(isequal(p1,p2) || isequal(p1,q2) || isequal(q1,p2) || isequal(q1,q2))
			return false
	end
    if ((o1 != o2) && (o3 != o4))
        return true
	end
    if (o1 == 0 && onSegment(p1, p2, q1))
		return true
	end
  
    if (o2 == 0 && onSegment(p1, q2, q1))
		return true
	end
  
    if (o3 == 0 && onSegment(p2, p1, q2))
		return true
	end
  
    if (o4 == 0 && onSegment(p2, q1, q2))
		return true; 
	end
	return false;
end

# ╔═╡ 1f3f2ac8-68d6-11eb-0ae6-972dd9115a49
intersects(Point(0,0),Point(5,5),Point(1,2),Point(1,5))

# ╔═╡ 726ab2b6-68d2-11eb-1077-0fd606de2d8e
intersects(Point(0,3.2421),Point(5.132,5.12),Point(1.132,2.123),Point(1.232,5.21))

# ╔═╡ 99167488-68c6-11eb-13db-f9cedf855276
begin
	deleteat!([1,3,1], [false,true,false])
end

# ╔═╡ b4113c3a-68c8-11eb-0982-33a1ce621b57
deleteat!([1,2,3],[false,true,false])

# ╔═╡ 754a0390-68ca-11eb-2c08-056dbccd852e
function inter_point(p1,p2,x1,x2)

	L1 = [ p1.y - p2.y, p2.x-p1.x, -(p1.x*p2.y - p2.x*p1.y) ]
	L2 = [ x1.y - x2.y, x2.x-x1.x, -(x1.x*x2.y - x2.x*x1.y) ]
    D  = L1[1] * L2[2] - L1[2] * L2[1]
    Dx = L1[3] * L2[2] - L1[2] * L2[3]
    Dy = L1[1] * L2[3] - L1[3] * L2[1]
    if(D != 0)
        x = Dx / D
        y = Dy / D
        return Point(x,y)
	end
end

# ╔═╡ d45665f6-68d7-11eb-2156-1f4a6b81bba1


# ╔═╡ aaa7cfd4-68cc-11eb-2f52-c9bf76acd4df
begin
aass = next_point.(genk.to,genk.angle,genk.isTip,genk.Si)
new_point= Point(1,2)
	
aass[1][2].x = new_point.x
end

# ╔═╡ 1b972c50-68d0-11eb-07b0-21c499061810
inter_point(Point(1,0),Point(1,5),Point(0,3),Point(5,3))

# ╔═╡ d94caad6-68df-11eb-28f6-c1df93195446
function branch_angle()
	return rand(Normal(1.39626, 0.174533))
end

# ╔═╡ d2bd7028-68f4-11eb-3fa0-a77fe1160848
function getI(c::CartesianIndex{2})
	return c[2]
end

# ╔═╡ 6e220a2a-68f4-11eb-2a32-17b9d6b1fef2
begin
ar = [-1 -1 1 1 -1]
getI.(findall(x -> x == 1,ar))
end

# ╔═╡ d41d5e40-68f6-11eb-0274-1ddf21a46161
function branch(p::Point, Si::Float64)
	if Si < 10^(-11)
		return false
	else
		if(bernoulli(Si*t*10^3))
			return true
		end
	end
	return false
end

# ╔═╡ eddde038-68f7-11eb-3368-6760f8a0763b
function shares_endpoint(p1,p2,x1,x2)
	if ((isequal(p1,x1) && isequal(p2,x2)) || (isequal(p1,x2) && isequal(p2,x1)))
			return false
	end
	return (isequal(p1,x1)
		|| isequal(p1,x2) 
		|| isequal(p2,x1) 
		|| isequal(p2,x2)) 
	
end

# ╔═╡ df7fd6dc-683d-11eb-2a9a-9feba787dd47
function step!(fungi::Fungi)
	# TIP EXTENSION
	s = next_point.(fungi.to,fungi.angle,fungi.isTip,fungi.Si)
	fungi.isTip = fungi.isTip .* 0
	s = filter(x -> x != -1, s)
	indices = []
	
	try
	append!(indices, getI.(findall(x -> x==1,s)))
	catch
		@warn "blabla"
	end
	s = filter(x -> x != 1, s)
	fungi.isTip[indices] .= 1
	
	
	# BRANCHING
#	for i in 1:length(fungi.from)
#		if(branch(fungi.from[i], fungi.Si[i]))
#			ϕ = branch_angle()
#			p = fungi.from[i]
#			
#			if length(s) == 0
#				s= [(p, Point(p.x + cos(ϕ), p.y +  sin(ϕ)), ϕ)]
#			else	
#				push!(s, (p, Point(p.x + cos(ϕ), p.y +  sin(ϕ)), ϕ))
#			end
#			
#		end
#	end
	
	# TRANSLOCATION
	dif = zeros(size(fungi.from))
	for i in 1:length(fungi.from)
		for j in 1:(length(fungi.from))
			if(shares_endpoint(fungi.from[i],fungi.to[i],fungi.from[j],fungi.to[j]))
				dif[i] +=  (fungi.Si[j] - fungi.Si[i])/(1e-4)
			end
		end
	end
	dif.*=t
	fungi.Si += dif
	## change the substrate!
	
	# UPTAKE
	for i in 1:length(fungi.from)
		xcoord = Int( floor(getx(fungi.to[i])/1e-4))
		ycoord = Int(floor(gety(fungi.to[i])/1e-4))
		ΔSi = 600*t*(fungi.Si[i]/(4*1e-9+fungi.Si[i]))
			 # Se[xcoord, ycoord ]
		#if ΔSi > Se[xcoord, ycoord]
			#fungi.Si[i] += Se[xcoord, ycoord]
			#Se[xcoord, ycoord] = 0
		#else
		#	fungi.Si[i] += ΔSi
			#Se[xcoord, ycoord] -= ΔSi
		#end
	end
	
	
	
	## INTERSECTION
	new_tips = ones(length(s))
	
	for i in 1:length(s)
		for j in 1:(length(fungi.to))
			p1 = s[i][1]
			p2 = s[i][2]
			x1 = fungi.from[j]
			x2 = fungi.to[j]
			if(intersects(p1,p2,x1,x2))
				new_p = inter_point(p1,p2,x1,x2)
				s[i][2].x = new_p.x
				s[i][2].y = new_p.y
				new_tips[i]=0
				break
			end
		end
	end

	append!(fungi.from, getFirst.(s))
	append!(fungi.to, getSecond.(s))
	append!(fungi.angle, getThird.(s))
	append!(fungi.isTip,new_tips )
	append!(fungi.Si, zeros(Float64,size(s)))
	
	return fungi
end

# ╔═╡ 96846800-6845-11eb-0a66-a7fb85b055bf
function getpoints(hypha::Hypha)
	hypha.points
end

# ╔═╡ 8d5f2f56-68c4-11eb-2a5c-2572d0e33ff4
genk

# ╔═╡ 400c487e-68d8-11eb-10f8-f719cf876699
intersects(Point(48.0,51.0),Point(47.13,51.5),Point(49.0,51.0),Point(48.0,51.00))

# ╔═╡ 791b8654-68ea-11eb-3483-0bf062d2508b
scatter([1],[1])


# ╔═╡ 01022fa4-68bf-11eb-10c4-c5c8049ff714
simulation

# ╔═╡ eaefdc8c-68a7-11eb-1cd5-575c2445e9ce
time

# ╔═╡ 730868de-6845-11eb-02b7-3b3a8f9f692a
function plot_fungi(fungi::Fungi)
	rng = range(HSV(0,1,1), stop=HSV(-360,1,1), length=200)
	sc = plot(50,50, legend=false)
	for i in 1:length(fungi.from)
		plot!([getx(fungi.from[i]), getx(fungi.to[i])], [gety(fungi.from[i]), gety(fungi.to[i])])
	end
	return sc
end

# ╔═╡ 77f4c4ba-6843-11eb-3920-636754a1d827
begin
	step!(genk)
	##@bind time Slider(1:length(getpoints(genk.hyphae[1])),default=1)
	plot_fungi(genk)
end


# ╔═╡ d9d7474c-68d3-11eb-292a-c5a466ab3e4b
begin
	plot_fungi(genk)
end

# ╔═╡ dea5193e-68c4-11eb-2a88-2b2c0ffca2e2
plot_fungi(simulation[i])

# ╔═╡ f3388864-68aa-11eb-3682-75dc9b9f5f4f
@gif for t in 1:length(getpoints(genk.hyphae[1]))
	plot_fungi(genk,t)
end

# ╔═╡ f0fd1926-68b3-11eb-3c50-1defc8f0b383
plot_fungi(genk)

# ╔═╡ 58cd588a-6851-11eb-1451-5bf9eb967c5f
md"
# Branching
Branching angles have been shownto be normally distributed in fungal mycelia (Hutchinson et al., 1980) and so angles are drawnfrom a Normal distributionN(μ, σ2).
"

# ╔═╡ Cell order:
# ╟─ec133ca2-6763-11eb-002e-45a5d78c3605
# ╟─c9c34aa6-675f-11eb-1fc5-196e34a58196
# ╟─b6da61c0-6775-11eb-1352-cf2f82d47fea
# ╟─0e794288-676d-11eb-0d34-1fc7c09880d5
# ╟─49420896-6761-11eb-025d-c7b05f2e0f08
# ╟─f6f8dca2-6764-11eb-1ea4-11949b018219
# ╟─1dae4998-6765-11eb-3906-5158627996f5
# ╟─a2efe0d2-6765-11eb-2249-c50e7ca75cd1
# ╟─62ff5ea8-6765-11eb-1ae9-e597834764df
# ╟─c4826416-6767-11eb-0e61-e912d529bebc
# ╟─19726662-6761-11eb-1d7a-df42ac99ef1c
# ╟─54503b68-6768-11eb-369c-9d1714c8a21c
# ╟─70cbb1c8-6768-11eb-36c4-b7f1e44ceb9d
# ╠═67f36fa8-676d-11eb-3055-bd1f93a5c936
# ╟─452e74e2-6772-11eb-347f-3556436653ca
# ╟─7a738940-6776-11eb-39e7-c33e9301faf3
# ╟─a1719b38-6774-11eb-01fd-511391d46d16
# ╟─aadd29bc-6774-11eb-02c6-59f021e86c6a
# ╟─35dfd5e6-6775-11eb-23d6-052efdcd1dd1
# ╟─8cea6466-6836-11eb-2704-797ca64396a1
# ╟─c7ecc5ca-6775-11eb-0da6-09045374100c
# ╟─7af43442-677a-11eb-2a4f-71ca69516b35
# ╟─e4efbe66-677a-11eb-105b-819ebbd28990
# ╠═4c2ed13a-6851-11eb-3fcb-c79d02727c89
# ╠═78c91270-6827-11eb-250a-6f5be3f36820
# ╠═ec175d6c-6828-11eb-2f6f-c75c5b18c3ef
# ╠═6892328c-68ff-11eb-2a85-4f5aa9d53b32
# ╠═ed19a464-68ea-11eb-1d13-cd2480f4186c
# ╠═1ea0c1a8-683d-11eb-0efa-89c527035050
# ╠═1c533af0-6829-11eb-2f2d-ab894ae524ca
# ╠═4d04ef40-68b0-11eb-2781-c9d566b2ece5
# ╠═55f94e50-6826-11eb-2a77-713890850b60
# ╠═2674f916-68d2-11eb-0066-1fcc4722e530
# ╠═d1293282-682f-11eb-1e05-71185e0dcc8b
# ╟─df20ac7a-6830-11eb-077f-5d0b553551c2
# ╟─f6aea732-6830-11eb-1fbe-bfd1f16815a3
# ╠═01e89ad2-684e-11eb-1a44-fde170bc69d3
# ╠═8babf7d4-684e-11eb-1cc5-8d13e27b32b3
# ╠═b9127852-68f5-11eb-0b11-172f8bb57677
# ╠═7cc884ac-68fc-11eb-0500-b3e0622146dc
# ╠═469a53c2-68fe-11eb-1064-8d806c382970
# ╟─7651344a-68f5-11eb-17be-d7573461d1a8
# ╠═ac8dd314-6847-11eb-3f2e-1163d16252d6
# ╠═398fdf3a-6854-11eb-3041-b9130d132483
# ╠═e9e11d3c-6853-11eb-0b2f-c310a84e7f9e
# ╠═f3df1314-68b4-11eb-2d92-6b515469e676
# ╟─d5796fb4-68b9-11eb-021a-e537d93880d7
# ╟─f6034e94-68b9-11eb-3586-d18d276bd9be
# ╟─331a2554-68bb-11eb-0fb9-99fc2c426c1a
# ╠═c097e286-68b9-11eb-2799-333f9b0cd757
# ╠═d6eac3f2-68be-11eb-026d-6729ab9283b1
# ╠═1c596458-68c8-11eb-3993-29610d466d9a
# ╠═00acdf7c-68d3-11eb-273a-61da7ca2de6f
# ╠═6d712a94-68c6-11eb-2881-1b394bde84a4
# ╠═1f3f2ac8-68d6-11eb-0ae6-972dd9115a49
# ╠═726ab2b6-68d2-11eb-1077-0fd606de2d8e
# ╠═99167488-68c6-11eb-13db-f9cedf855276
# ╠═b4113c3a-68c8-11eb-0982-33a1ce621b57
# ╠═754a0390-68ca-11eb-2c08-056dbccd852e
# ╠═d45665f6-68d7-11eb-2156-1f4a6b81bba1
# ╠═aaa7cfd4-68cc-11eb-2f52-c9bf76acd4df
# ╠═1b972c50-68d0-11eb-07b0-21c499061810
# ╠═d94caad6-68df-11eb-28f6-c1df93195446
# ╠═d2bd7028-68f4-11eb-3fa0-a77fe1160848
# ╠═6e220a2a-68f4-11eb-2a32-17b9d6b1fef2
# ╠═d41d5e40-68f6-11eb-0274-1ddf21a46161
# ╠═eddde038-68f7-11eb-3368-6760f8a0763b
# ╠═df7fd6dc-683d-11eb-2a9a-9feba787dd47
# ╟─96846800-6845-11eb-0a66-a7fb85b055bf
# ╟─3a5fc9dc-68ab-11eb-3a40-9fc97f9b2ebe
# ╠═8d5f2f56-68c4-11eb-2a5c-2572d0e33ff4
# ╠═400c487e-68d8-11eb-10f8-f719cf876699
# ╠═77f4c4ba-6843-11eb-3920-636754a1d827
# ╠═d9d7474c-68d3-11eb-292a-c5a466ab3e4b
# ╠═791b8654-68ea-11eb-3483-0bf062d2508b
# ╠═01022fa4-68bf-11eb-10c4-c5c8049ff714
# ╠═eaefdc8c-68a7-11eb-1cd5-575c2445e9ce
# ╠═dea5193e-68c4-11eb-2a88-2b2c0ffca2e2
# ╠═f3388864-68aa-11eb-3682-75dc9b9f5f4f
# ╟─730868de-6845-11eb-02b7-3b3a8f9f692a
# ╠═f0fd1926-68b3-11eb-3c50-1defc8f0b383
# ╠═58cd588a-6851-11eb-1451-5bf9eb967c5f
