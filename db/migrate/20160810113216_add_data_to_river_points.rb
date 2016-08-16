class AddDataToRiverPoints < ActiveRecord::Migration
  def change
    add_column :river_points, :river, :string
    add_column :river_points, :branch, :integer
    add_column :river_points, :point, :integer
    add_column :river_points, :latitude, :float, :limit => 30, :scale => 12
    add_column :river_points, :longitude, :float, :limit => 30, :scale => 12
  end
end
