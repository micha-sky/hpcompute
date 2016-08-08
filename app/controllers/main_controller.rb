class MainController < ApplicationController
  require 'csv'
  def index

  end

  def map

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

  def read_outlocs
    filename = 'outlocs/outloc.0001'

    CSV.foreach(filename, { :headers => true, :col_sep => "\t", :skip_blanks => true }) do |row|
      puts row
    end

    render :text => File.read(filename)

  end

  private
  def cache_results
    if File.exists?('pcolour.png')
      true
    end
  end


end
