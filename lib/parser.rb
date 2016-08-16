module Parser
  require 'csv'
  extend self

  # @return [hash]
  def read_outlocs(filename, river, branch, point)
    result = {}
    times = []
    thalwegs = []
    depths = []
    discharges = []
    file = CSV.open(filename, { :headers => true, :col_sep => "\t", :skip_blanks => true })
    file.drop(3).each do |row|
      times << row[0].split(' ')[0]
      thalwegs << row[0].split(' ')[3]
      depths << row[0].split(' ')[5]
      discharges << row[0].split(' ')[6]
    end

    result[:river] = river
    result[:branch] = branch
    result[:point] = point
    result[:times] = times
    result[:discharges] = discharges
    result[:thalwegs] = thalwegs
    result[:depths] = depths

    result

  end
end