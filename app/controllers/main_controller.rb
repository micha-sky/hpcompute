class MainController < ApplicationController
  def index

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

  private
  def cache_results
    if File.exists?('pcolour.png')
      true
    end
  end

  def read_outlocs

  end
end
