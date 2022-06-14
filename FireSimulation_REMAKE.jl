using GLMakie
using Random, Distributions
using LinearAlgebra
using Combinatorics

const N = 6000
const Window = 100
const Res = 800

const Size = 0.40
const max_size = 0.6
const min_size = 0.35


const Tmax = 1000
const dt = 0.02


const spawn_rate = 100
const spawn_nb = 0

const Wall = true
const Roof = true

const size_coef = 10.84 #5.42*2

const sigma = 0.04
const Thermostat_temp = 300
const Thermostat_time = 0
const temp_threshold = 250
const temp_speed_threshold = 250
const temp_force = 150
const temp_speed_force = 0.8

const heat_radiation = 0.05

const compression = 0.5



function init()
    #INIT ARRAYS ============================================
    POS_LIST = Observable(GLMakie.Point2f0[])
    OLD_POS_LIST = Observable(GLMakie.Point2f0[])

    PRESSURE_LIST = Observable(Float64[])
    PRESSURE_MEAN = Observable(Float64[])
    TEMP_LIST = Observable(Float64[])
    #display
    SPRITE_LIST = Observable(Float64[])
    COLOR_LIST = Observable(Float64[])


    FPS_title = Observable("FPS : XX")
    for _ in 1:N
        x,y = Window*rand(), Window*rand()*compression
        push!(POS_LIST[], Point2f0(x, y))
        push!(OLD_POS_LIST[], Point2f0(x, y)) 
        push!(SPRITE_LIST[], Size*size_coef)
        push!(PRESSURE_LIST[], 0)
        push!(TEMP_LIST[], 0)
        push!(COLOR_LIST[],rand())
    end 
    #INIT PLOT ==============================================
    fig = Figure()

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()

    axtop = Axis(ga[1, 1], title = FPS_title)
    axbot = Axis(gb[1, 1])
    rowsize!(fig.layout, 1, Auto(5))
    #poly!(axtop,[Rect(0,0,Window,Window) for i in 1:2], color = 2:-1:1, colormap = :thermal)
    scatter!(axtop,POS_LIST, markersize = SPRITE_LIST ,color = COLOR_LIST,colormap = :thermal)
    

    limits!(axtop,0,Window,0,Window)
    axtop.xticks = 0:5:Window
    axtop.yticks = 0:5:Window

    axtop.aspect = 1
    axbot.aspect = 5
    lines!(axbot,PRESSURE_MEAN)
    limits!(axbot,0,Tmax,0,3) #Pressure
    screen = display(fig)
    
    resize!(screen, Res, Res)
    println("Init finished")
    SIMU(fig,POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST,PRESSURE_MEAN,COLOR_LIST,FPS_title)
end


function SIMU(fig,POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST,PRESSURE_MEAN,COLOR_LIST,FPS_title)
    t1,t2 = 0,1
    Thermostat = false
    PRESSURE_LIST_tempo = Float64[]
    GLMakie.record(fig, "FireSimu_REMAKE.mp4", 1:Tmax) do t
        TotalTime = t2 - t1
        t1 = time()
        FPS_title[] = "FPS : " * string(round(1/TotalTime,digits = 2))
        PRESSURE_LIST_tempo = zeros(length(POS_LIST[]))

        if t > Thermostat_time 
            Thermostat = true
        end

        ComputeTime = @elapsed POS_LIST, OLD_POS_LIST,PRESSURE_LIST_tempo = update(POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST_tempo,dt,Thermostat)
        RenderTime = TotalTime - ComputeTime
        if (t%spawn_rate==0 && t!=0)
            spawn(spawn_nb,POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST_tempo)
        end
        lambda = 0.05
        a = (max_size-min_size)/300
        
          
        COLOR_LIST[] = 300 ./(1 .+exp.(-lambda*(TEMP_LIST[] .- 200)))
        SPRITE_LIST[]=  (a.*TEMP_LIST[].+min_size)*size_coef 


        PRESSURE_LIST[] = PRESSURE_LIST_tempo
        notify.((POS_LIST,PRESSURE_LIST,SPRITE_LIST, COLOR_LIST, PRESSURE_MEAN))
        push!(PRESSURE_MEAN[],mean(PRESSURE_LIST[]))
        
        println(t,"- N: ", length(POS_LIST[]), "- ComputeTime : ", round(ComputeTime*1000,digits=2)," ms", "| RenderTime : ", round(RenderTime*1000,digits=2)," ms -- Pressure: ", round(mean(PRESSURE_LIST[]),digits=2))
        t2 = time()
        
    end
    return POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST,COLOR_LIST
end



function update(POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST,dt,Thermostat)

    POS_LIST,OLD_POS_LIST = updatePositions(POS_LIST,OLD_POS_LIST,TEMP_LIST,dt)
    POS_LIST,OLD_POS_LIST,TEMP_LIST = applyConstraint(POS_LIST,OLD_POS_LIST, TEMP_LIST,Thermostat)
    POS_LIST,PRESSURE_LIST, TEMP_LIST = applyCollision(POS_LIST,PRESSURE_LIST,TEMP_LIST,Thermostat)
    
    return POS_LIST,OLD_POS_LIST,PRESSURE_LIST
end

function updatePositions(POS_LIST,OLD_POS_LIST,TEMP_LIST,dt)
    for i in 1:length(POS_LIST[])
        pos,old_pos = POS_LIST[][i],OLD_POS_LIST[][i]
        temp = TEMP_LIST[][i]
        POS_LIST[][i],OLD_POS_LIST[][i] = updatePos(pos,old_pos,temp,dt)
    end
    return POS_LIST,OLD_POS_LIST
end

function updatePos(pos,old_pos,temp,dt)
    vel = pos - old_pos
    old_pos = pos
    acc = getAcc(temp)
    pos = pos + vel + acc*dt*dt
    return pos,old_pos
end

function getAcc(temp)
    temp_acc = [0,0]
    if temp > temp_threshold
        temp_acc = [0,temp_force]
    end
    return [0,-10.0] + temp_acc
end

function applyCollision(POS_LIST,PRESSURE_LIST,TEMP_LIST,Thermostat)

    for i in 1:length(POS_LIST[])

        pos1 = POS_LIST[][i]
        for j in 1:length(POS_LIST[])
            if (i==j)
                continue
            end
            pos2 = POS_LIST[][j]
            axis = pos1 - pos2
            dist = norm(axis)
            if (dist < 2*Size && dist > 0.001)
                u = axis/dist
                delta = 2*Size - dist
                coef = 0.5
                POS_LIST[][i] += 0.5*delta*u*coef #Put viscosity back
                POS_LIST[][j] -= 0.5*delta*u*coef
                PRESSURE_LIST[i] += delta
                PRESSURE_LIST[j] += delta

                if Thermostat
                    
                    TEMP_LIST[][i],TEMP_LIST[][j] = heat_transfer(TEMP_LIST[][i],TEMP_LIST[][j], sigma)
                end
            end
        end
    end
    return POS_LIST,PRESSURE_LIST,TEMP_LIST

end
function applyConstraint(POS_LIST, OLD_POS_LIST,TEMP_LIST,Thermostat)

    for i in 1:length(POS_LIST[])
        new_pos = POS_LIST[][i]
        radius = Size

        TEMP_LIST[][i] = max(0,TEMP_LIST[][i] - heat_radiation)

        if new_pos[2]  - radius < 0 #bouce on ground
            new_pos = Point2f0(new_pos[1],radius)
            if Thermostat
                temp = TEMP_LIST[][i]
                TEMP_LIST[][i], _ = heat_transfer(temp,Thermostat_temp, 0.25)
                if TEMP_LIST[][i] > temp_speed_threshold
                    OLD_POS_LIST[][i] = Point2f0(new_pos[1], new_pos[2] - temp_speed_force)
                end
            end
        elseif new_pos[2] + radius > Window && Roof #ROOF
            new_pos = Point2f0(new_pos[1] ,Window - radius)
        end
        if new_pos[1] -radius < 0 && Wall#LEFT
            new_pos = Point2f0(radius,new_pos[2])

        elseif new_pos[1] + radius> Window && Wall#RIGHT
            new_pos = Point2f0(Window - radius,new_pos[2])
        end
        POS_LIST[][i] = new_pos
    end
    
    return POS_LIST,OLD_POS_LIST,TEMP_LIST
end

function heat_transfer(temp1,temp2, sigma)
    m = (temp1-temp2)/2
    temp1 -= m*sigma
    temp2 += m*sigma
    return temp1,temp2
end

function get_size(temp)
    a = 300/(max_size-min_size)
    return a*temp+min
    
end


function spawn(n,POS_LIST,OLD_POS_LIST,SPRITE_LIST,TEMP_LIST,PRESSURE_LIST)
    println("spawn")
    for _ in 1:n
        x,y = 50.0,75.0 #rand()*Window*1 #rand()*window, rand()*window
        v = [0.0,-0.0]
        push!(POS_LIST[], Point2f0(x,y))
        push!(OLD_POS_LIST[], Point2f0(x,y) - v)
        push!(PRESSURE_LIST[], 0)
        push!(SPRITE_LIST[],Size*size_coef)
        push!(TEMP_LIST[],0)
    end
    return POS_LIST,OLD_POS_LIST,PRESSURE_LIST
end
init()