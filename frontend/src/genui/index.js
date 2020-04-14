import RoutedPage from "./components/RoutedPage";
import ResponsiveGrid from "./components/ResponsiveGrid";
import TabWidget from './components/TabWidget';
import TaskAwareComponent from './components/tasks/TaskAwareComponent';
import TaskBadgeGroup from './components/tasks/TaskBadgeGroup';
import TaskProgressBar from './components/tasks/TaskProgressBar';
import ComponentWithObjects from './components/ComponentWithObjects';
import FieldErrorMessage from './components/forms/FieldErrorMessage';
import ProjectItemSubTitle from './components/ProjectItemCardSubtitle';
import {TableHeaderFromItems, TableDataFromItems} from './components/tables/Tables';
import DownloadFile from './components/DownloadFile';
import LiveObject from './components/LiveObject';
import ComponentWithResources from './components/ComponentWithResources';
import ModelGrid from './components/models/ModelGrid';
import ModelsPage from './components/models/ModelsPage';
import ModelCardNew from './components/models/ModelCardNew';
import ModelCard from './components/models/ModelCard';
import ModelInfoTab from './components/models/tabs/ModelInfo';
import ModelPerformanceTab from './components/models/tabs/ModelPerf';
import GenericMolSetCard from './components/compounds/GenericMolSetCard';
import GenericMolSetGrid from './components/compounds/GenericMolSetGrid';
import GenericNewMolSetCard from './components/compounds/GenericNewMolSetCard';
import CompoundsPage from './components/compounds/CompoundsPage';
import MolSetTasks from './components/compounds/MolSetTasks';
import FileUpload from './components/forms/FileUpload';
import FormikModelUploadForm from './components/models/FormikModelUploadForm';
import {groupByMolset, filterProviders, resolve, smoothScrollToTop, IDsToResources, groupBy, scrollTo} from './utils'
import {MoleculeImage, MoleculePic} from './components/compounds/details/MoleculeImage'
import {MoleculeMetadata, DataPair} from './components/compounds/details/MoleculeMetadata';
import ComponentWithPagedResources from './components/ComponentWithPagedResources';
import MoleculeActivityProvider from './components/compounds/details/MoleculeActivityProvider';
import {ActivitiesByTypeTabView, ActivitySetTabView, ActivitiesTable, ActivitySetFlatView, ActivitiesByTypeFlatView} from './components/compounds/details/ActivityViews'
import TabWidgetSmart from './components/TabWidgetSmart';
import MolsToMolSetGroups from './components/compounds/summaries/MolsToMolSetGroups';
import MolSetsTabs from './components/compounds/summaries/MolSetsTabs';
import CompoundList from './components/compounds/summaries/CompoundList';
import ApiResourcePaginator from './components/ApiResourcePaginator';
import CompoundListFromAPI from './components/compounds/summaries/CompoundListFromAPI';
import MolsInMolSetList from './components/compounds/tabs/MolsInMolSetList';
import ModelPreds from './components/models/tabs/predictions/ModelPreds';
import GroupedViolinPlot from './components/GroupedViolinPlot';
import ActivitiesAggregator from './components/compounds/summaries/ActivitiesAggregator';
import CompoundOverview from './components/compounds/summaries/CompoundOverview';
import MoleculePropsProvider from './components/compounds/details/MoleculePropsProvider';
import {PropertiesTable} from './components/compounds/details/PropertyViews'
import SimpleDropDownToggle from './components/SimpleDropDownToggle';
import {ActivitySetStatsTable, MolsetActivitiesSummary} from './components/compounds/summaries/ActivitySetStatsTable';

// TODO: structure this list in a more sensible way
export {
  RoutedPage,
  ResponsiveGrid,
  TabWidget,
  TabWidgetSmart,
  TaskAwareComponent,
  TaskBadgeGroup,
  TaskProgressBar,
  ComponentWithObjects,
  FieldErrorMessage,
  ProjectItemSubTitle,
  TableHeaderFromItems,
  TableDataFromItems,
  resolve,
  DownloadFile,
  LiveObject,
  ComponentWithResources,
  ModelGrid,
  ModelsPage,
  ModelCardNew,
  ModelCard,
  ModelInfoTab,
  ModelPerformanceTab,
  GenericMolSetCard,
  GenericMolSetGrid,
  GenericNewMolSetCard,
  CompoundsPage,
  MolSetTasks,
  FileUpload,
  FormikModelUploadForm,
  groupByMolset,
  MoleculeImage,
  MoleculePic,
  MoleculeMetadata,
  DataPair,
  filterProviders,
  ComponentWithPagedResources,
  MoleculeActivityProvider,
  ActivitiesByTypeTabView,
  ActivitySetTabView,
  ActivitiesTable,
  MolsToMolSetGroups,
  MolSetsTabs,
  CompoundList,
  smoothScrollToTop,
  ApiResourcePaginator,
  IDsToResources,
  CompoundListFromAPI,
  MolsInMolSetList,
  ModelPreds,
  groupBy,
  GroupedViolinPlot,
  ActivitiesAggregator,
  CompoundOverview,
  MoleculePropsProvider,
  PropertiesTable,
  ActivitySetFlatView,
  SimpleDropDownToggle,
  ActivitiesByTypeFlatView,
  ActivitySetStatsTable,
  MolsetActivitiesSummary,
  scrollTo
}