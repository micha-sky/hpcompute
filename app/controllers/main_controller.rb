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
    head 400 if params.blank?

    river = params[:river]
    tmp_dir = "tmp/#{river}-#{timestamp}"

    Dir.mkdir(tmp_dir)
    river_files = %W(models/rivtoxClean/rivtox_new.exe models/rivtoxClean/input_#{river} models/rivtoxClean/extraction_points_#{river}.xy)
    river_files.each do |file|
      FileUtils.cp(file, tmp_dir)
    end

    render :text => tmp_dir

    Dir.chdir(tmp_dir) do
      File.rename("input_#{river}", 'input')
      command = `./rivtox_new.exe`
      command
    end

  end

  def get_results
    river = params[:river]
    branch = params[:branch]
    point = params[:point]
    path = params[:path]

    filename = File.join(path, 'outlocs/outloc.0001')

    result = Parser.read_outlocs(filename, river, branch, point)

    render :json => result.to_json

  end

  private
    def timestamp
      Time.now.to_i
    end
end
