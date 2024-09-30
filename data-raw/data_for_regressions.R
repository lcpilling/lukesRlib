
# example data from the excellent folk at https://github.com/ErikinBC/SurvSet
example_data = readr::read_csv("https://github.com/ErikinBC/SurvSet/raw/main/SurvSet/_datagen/output/Framingham.csv")
#example_data = rbind(example_data,example_data)
example_data = example_data |>
  dplyr::select(-pid) |>
  dplyr::rename(
         sbp=num_sbp,
         dbp=num_dbp,
         scl=num_scl,
         age=num_age,
         bmi=num_bmi,
         sex=fac_sex,
         month=fac_month)

example_data = example_data |> dplyr::mutate(
  
  bmi_cat = dplyr::case_when(
    bmi <  25 ~ 0,
    bmi >= 25 & bmi < 30 ~ 1,
    bmi >= 30 ~ 2),
  bmi_cat = haven::labelled(bmi_cat, c(Normal = 0, Overweight = 1, Obese = 2), label="BMI (categories)"),
  
  sex = as.numeric(dplyr::if_else(
    sex == "M", 1, 0
  )),,
  sex = haven::labelled(sex, c(Female = 0, Male = 1), label="Sex")
  
)


# test
#lukesRlib::tidy_ci( glm( sbp ~ age + as.factor(sex) + as.factor(bmi_cat) + scl, data = example_data ) )
#lukesRlib::tidy_ci( glm( sex ~ age + as.factor(bmi_cat) + sbp + scl, data = example_data, family = binomial ) )
#lukesRlib::tidy_ci( coxph( Surv(time, event) ~ age + as.factor(sex) + as.factor(bmi_cat) + sbp + scl, data = example_data ) )
#lukesRlib::get_assoc( y="event", x=c("bmi","sbp","scl"), z="age+as.factor(sex)", d=example_data, model="logistic" )
#lukesRlib::get_assoc( y="Surv(time, event)", x=c("bmi","sbp","scl"), z="age+as.factor(sex)", d=example_data, model="coxph" )


# save
usethis::use_data(example_data, overwrite = TRUE, compress = 'xz')

