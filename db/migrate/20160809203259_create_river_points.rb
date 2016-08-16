class CreateRiverPoints < ActiveRecord::Migration
  def change
    create_table :river_points do |t|

      t.timestamps null: false
    end
  end
end
