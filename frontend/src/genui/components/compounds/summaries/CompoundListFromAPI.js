import { IDsToResources } from '../../../utils';
import { ComponentWithResources, CompoundList } from '../../../index';
import ApiResourcePaginator from '../../ApiResourcePaginator';
import React from 'react';

export default function CompoundListFromAPI(props) {
  const activitySetsIDs = props.activitySetsIDs ? props.activitySetsIDs : [];

  const url = new URL(`${props.molset.id}/molecules/`, props.apiUrls.compoundSetsRoot);
  if (activitySetsIDs.length > 0) {
    const resourcesDef = IDsToResources(props.apiUrls.activitySetsRoot, activitySetsIDs);
    return (
      <ComponentWithResources
        {...props}
        definition={resourcesDef}
      >
        {
          (allLoaded, activitySets) => {
            return allLoaded ? (
              <ApiResourcePaginator
                {...props}
                url={url}
              >
                {
                  currentItems => {
                    return (
                      <CompoundList
                        {...props}
                        activitySets={activitySets}
                        mols={currentItems}
                        paginate={false}
                      />
                    )
                  }
                }
              </ApiResourcePaginator>
            ) : <div>Loading...</div>
          }
        }
      </ComponentWithResources>
    )
  } else {
     return (
       <ApiResourcePaginator
         {...props}
         url={url}
       >
         {
           currentItems => {
             return (
               <CompoundList
                 {...props}
                 activitySets={[]}
                 showInfo={true}
                 showActivities={false}
                 mols={currentItems}
                 paginate={false}
               />
             )
           }
         }
       </ApiResourcePaginator>
     )
  }
}