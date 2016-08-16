class MainController < ApplicationController
  def index

  end

  def map

  end

  def get_points
    river = params[:river]
    points = RiverPoint.where(:river => river)

    render :json => points.to_json
  end


  def execute
    command = `make plots -C models/jacobi`
    render :text => command.to_s
  end

  def executeRivtox
    command = `./models/Kanev_new/a.out`
    render :text => command.to_s
  end

  def executeRivtox2
    command = `./models/rivtoxClean/rivtox_new.exe`
    render :text => command.to_s
  end

  def get_results
    river=params[:river]
    branch=params[:branch]
    point=params[:point]

    filename = 'outlocs/outloc.0001'

    result = Parser.read_outlocs(filename, river, branch, point)

    render :json => result.to_json

  end
end
