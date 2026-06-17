from django.urls import path, include
from rest_framework.routers import DefaultRouter
from . import api_views

# 创建路由器并注册我们的视图集
router = DefaultRouter()
router.register(r'hormones', api_views.HormoneViewSet)
router.register(r'literature', api_views.LiteratureViewSet)
router.register(r'analyses', api_views.AnalysisHistoryViewSet)

urlpatterns = [
    # API根路径，自动生成路由
    path('', include(router.urls)),
    
    # 基因分析API
    path('analyze/', api_views.GeneAnalysisAPIView.as_view(), name='gene-analysis'),
    
    # 用户API
    path('users/me/', api_views.CurrentUserView.as_view(), name='current-user'),
    path('users/register/', api_views.UserRegistrationAPIView.as_view(), name='user-register'),
    path('users/login/', api_views.UserLoginAPIView.as_view(), name='user-login'),
    
    # 异步任务API
    path('tasks/', api_views.TaskListView.as_view(), name='task-list'),
    path('tasks/<str:task_id>/', api_views.TaskDetailView.as_view(), name='task-detail'),
    
    # API文档
    path('docs/', include('rest_framework.urls', namespace='rest_framework'))
]    